#!/usr/local/bin/Rscript --vanilla

args = commandArgs(trailingOnly=TRUE)

sample              = args[1]  # name until ".barcode.cell.distribution.tsv"
file                = args[2] # path and file name "~/T23_8_lib_00011AAE_ACGGAACA-ACGAGAAC.barcode.cell.distribution.tsv" 
doublets            = args[3] # path and file name doblets
path                = args[4] # path to save the results
NonHhaI             = args[5] # path and file name control amplicons, csv with amplicons names
CpG.always.unmeth.B = args[6] # path and file name always digest control amplicons, csv with amplicons names
undigested          = args[7] # New argument: "YES", "yes", "Yes" or "NO"

# Extract filename prefix (before first dot)
name    = unlist(strsplit(basename(file), "\\."))[1]

# Ensure path ends with "/"
if (!grepl("/$", path)) {
  path <- paste0(path, "/")
}

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

# Read data
data <- read.table(file, header = TRUE, sep = "\t", row.names = 1)

# Create Seurat object
seurat <- CreateSeuratObject(t(data), 
                             project = "epimutations", 
                             assay = "DNAm.counts")

# Convert to binary data
binary_data  <- as.data.frame(ifelse(t(data) > 0, 1, 0))
binary_assay <- CreateAssayObject(counts = binary_data)

seurat[["DNAm"]] <- binary_assay

# Read doublet information
doublets <- read.csv(doublets, row.names = 1)

# Ordering
doublets <- doublets[doublets$Barcode %in% colnames(seurat),]

# Add metadata
seurat[["sample"]]           <- sample
seurat[["doublets"]]         <- doublets$DoubletDetectionLabel
seurat[["doublets.score"]]   <- doublets$DoubletDetectionScore

# Read control amplicons
NonHhaI <- read.csv(NonHhaI, header = FALSE)[,1]

# Compute QC metrics for NonHhaI
seurat$NonHhaI <- colSums(seurat[["DNAm"]]@counts[rownames(seurat[["DNAm"]]@counts) %in% NonHhaI, ])

# NonHhaI.cutoff
NonHhaI.cutoff = as.numeric(0.75 * length(NonHhaI))  

# Normalize `undigested` argument to lowercase
undigested <- tolower(undigested)

# Apply different filtering rules based on `undigested`
if (undigested == "yes") {
  
  # Read control amplicons
  CpG.unmeth <- read.csv(CpG.always.unmeth.B, header = FALSE)[,1]
  
  # Compute QC metrics for NonHhaI
  seurat$CpG.unmeth <- colSums(seurat[["DNAm"]]@counts[rownames(seurat[["DNAm"]]@counts) %in% CpG.unmeth, ])
  
  
  seurat.subset <- subset(seurat, 
                          subset = CpG.unmeth < 1 &
                            NonHhaI > NonHhaI.cutoff & 
                            doublets == 0 &
                            nCount_DNAm.counts < 400000 &
                            nCount_DNAm.counts > 8000)
} else {
  seurat.subset <- subset(seurat, 
                          subset = NonHhaI > NonHhaI.cutoff & 
                            doublets == 0 &
                            nCount_DNAm.counts < 400000 &
                            nCount_DNAm.counts > 8000)
}



# Seurat v4 vs. v5 compatibility
if (packageVersion("Seurat") >= "5.0.0") {
  counts <- GetAssayData(seurat.subset, assay = "DNAm.counts", layer = "counts")
} else {
  counts <- GetAssayData(seurat.subset, assay = "DNAm.counts", slot = "counts")
}

# Save filtered count matrix
write.table(t(as.matrix(counts)), paste0(path, name, ".barcode.cell.distribution_filtered.tsv"), 
            sep= "\t", row.names = TRUE, col.names = TRUE, quote= FALSE)

barcodes = row.names(t(as.matrix(counts)))

# Save filtered count matrix
write.table(barcodes, paste0(path, name, ".good.barcodes.tsv"), row.names = FALSE, col.names = FALSE, quote= FALSE)


# Generate histograms 

tryCatch({
  pdf(paste0(path, name, ".QC_plots.pdf"), onefile = TRUE)
  
  plot_list <- list(
    ggplot(seurat.subset@meta.data, aes(x = NonHhaI)) +
      geom_histogram(bins = 100, fill = "#D6F6FF") +
      theme_classic() + ggtitle("Distribution NoHhaI amplicons"),
    
    ggplot(seurat.subset@meta.data, aes(x = nFeature_DNAm)) +
      geom_histogram(bins = 100, fill = "#E6B5B9") +
      theme_classic() + ggtitle("Distribution nAmplicons"),
    
    ggplot(seurat.subset@meta.data, aes(x = nCount_DNAm.counts)) +
      geom_histogram(bins = 100, fill = "#FCAE7C") +
      theme_classic() + ggtitle("Distribution nCounts"),
    
    ggplot(seurat.subset@meta.data, aes(x = doublets.score)) +
      geom_histogram(bins = 100, fill = "#FFE699") +
      theme_classic() + ggtitle("Distribution doublet score")
  )
  
  # Add CpG.always.unmeth.B plot **only if undigested == "yes"**
  
  if (undigested == "yes") {
    plot_list <- c(plot_list, list(
      ggplot(seurat.subset@meta.data, aes(x = CpG.unmeth)) +
        geom_histogram(bins = 100, fill = "#E2CBF7") +
        theme_classic() + ggtitle("Distribution CpG always unmeth")
    ))
  }
  
  # Print all plots to PDF
  for (p in plot_list) {
    print(p)
  }
  
  dev.off()
}, error = function(e) {
  message("Warning: An error occurred while generating plots. The filtering process completed successfully!")
})
