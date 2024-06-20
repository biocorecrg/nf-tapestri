#!/usr/local/bin/Rscript --vanilla

#sample specific parameters

args = commandArgs(trailingOnly=TRUE)

sample              = args[1]  # name until ".barcode.cell.distribution.tsv"
file                = args[2]
doublets            = args[3]
path                = args[4]
NonHhaI             = args[5]
CpG.always.unmeth.B = args[6]


## Libraries -------------------------------------------------------------------
require(Seurat)
require(ggplot2)
require(tidyr)
require(SeuratObject)
require(ggpubr)
require(dplyr)
require(tidyverse)
require(reshape2)

# Read data & Create Seurat object ---------------------------------------------
data <- read.table(paste0(path, file), 
                         header = T, sep = "\t")
    
seurat <- data %>% 
      t() %>% 
      CreateSeuratObject(project = "epimutations",
                       min.cells = 0,
                       min.features = 0,
                       assay = "DNAm.counts") 
  
# create assay with binary data
binary_data  <- as.data.frame(ifelse(data > 0, 1, 0)) %>%
    t() %>%
    CreateAssayObject(min.cells = 0,
                      min.features = 0)
 
seurat[["DNAm"]] <- binary_data
  
# reading doublets information---------------------------------------------------
 
doublets <- read.csv(paste0(path, doublets), 
                         row.names="Barcode") 
### ensure order
doublets <- doublets[rownames(doublets) %in% colnames(seurat),]
all(rownames(doublets) == colnames(seurat))
  
## Adding information to Seurat object------------------------------------------
seurat[["sample"]]           <- c(rep(sample, length(colnames(seurat))))
seurat[["doublets"]]         <- doublets$DoubletDetectionLabel
seurat[["doublets.score"]]   <- doublets$DoubletDetectionScore

# reading control amplicons
NonHhaI <- read.csv(NonHhaI, header = F)
CpG.always.unmeth.B <- read.csv(CpG.always.unmeth.B, header = F)

# add QC varibles to seurat object
seurat@meta.data$CpG.always.unmeth.B = colSums(seurat[["DNAm"]]@counts[rownames(seurat@assays$DNAm) %in% CpG.always.unmeth.B$V1, ])
seurat@meta.data$NonHhaI = colSums(seurat[["DNAm"]]@counts[rownames(seurat@assays$DNAm) %in% NonHhaI$V1, ])

# First QC plots
Idents(seurat) = "sample"

pdf(paste0(path, sample, ".plot.before.QC.pdf"), onefile = TRUE)

hist(seurat$CpG.always.unmeth.B, breaks = 100, main = "Distribuion CpG always unmeth", col = "#E2CBF7")
hist(seurat$NonHhaI, breaks = 100, main = "Distribution NoHhaI amplicons", col = "#D6F6FF")
hist(seurat$nFeature_DNAm, breaks = 100, main = "Distribution nAmplicons", col = "#E6B5B9")
hist(seurat$nCount_DNAm.counts, breaks = 100, main = "Distribution nCounts", col = "#FCAE7C")
hist(seurat$doublets.score, breaks = 100, main = "Distribuion doublet score", col = "#FFE699")

dev.off()

# define undigested amplicons cut off
NonHhaI.cutoff = 0.75 * length(NonHhaI)

#SUBSET SEURAT OBEJECT based on quality
seurat.subset <- subset(seurat,
                        subset = CpG.always.unmeth.B < 1 &
                          NonHhaI > NonHhaI.cutoff & 
                          doublets == 0 &
                          nCount_DNAm.counts < 400000 &
                          nCount_DNAm.counts > 8000)

Idents(seurat.subset) = "sample"
pdf(paste0(path, sample, ".plot.after.QC.pdf"), onefile = TRUE)

hist(seurat.subset$CpG.always.unmeth.B, breaks = 100, main = "Distribuion CpG always unmeth", col = "#E2CBF7")
hist(seurat.subset$NonHhaI, breaks = 100, main = "Distribution NoHhaI amplicons", col = "#D6F6FF")
hist(seurat.subset$nFeature_DNAm, breaks = 100, main = "Distribution nAmplicons", col = "#E6B5B9")
hist(seurat.subset$nCount_DNAm.counts, breaks = 100, main = "Distribution nCounts", col = "#FCAE7C")
hist(seurat.subset$doublets.score, breaks = 100, main = "Distribuion doublet score", col = "#FFE699")
dev.off()

print(seurat.subset)

# extract filtered count matrix

counts <- GetAssayData(seurat.subset, 
                       assay = "DNAm.counts", 
                       layer = "counts")

# Save filtered matrix
write.table(t(as.matrix(counts)), 
            paste0(path, sample, ".barcode.cell.distribution_filtered.tsv"), 
            sep= "\t", 
            row.names = TRUE, 
            col.names = TRUE, 
            quote= FALSE)

