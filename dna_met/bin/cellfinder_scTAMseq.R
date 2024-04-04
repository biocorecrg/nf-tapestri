#!/usr/bin/Rscript --vanilla

require(argparse)
ap <- ArgumentParser()
ap$add_argument("-f", "--file", action="store", help="The barcode distribution file as output of the MissionBio pipeline (tsv)")
ap$add_argument("-a", "--ampli", action="store", help="A file containing amplicons for cell selection")
ap$add_argument("-o", "--output", action="store", help="The output file name")
ap$add_argument("-c", "--cutoff", action="store", type="double", default=0.7, help="Integer number specifiying the number of amplicons that need to fulfill the coverage threshold")
opt <- ap$parse_args()

# V1 tends to perform better for us
vers.2 <- FALSE
# cut off for amplicons detection
val <- as.numeric(opt$cutoff)

# modified cellfinder function from MISSION BIO
resemble.cellfinder <- function(d){
  cut <- ifelse(vers.2,ncol(d) * 8,ncol(d) * 10)
  reads.per.cell <- apply(d,1,sum)
  print(paste("Removing",sum(reads.per.cell<cut),"cells due to amplicon cutoff"))
  d <- d[reads.per.cell>=cut,]
  m.ean <- mean(apply(d,1,mean))
  thres <- ifelse(vers.2,0.2*m.ean,min(0.2*m.ean,10))
}

# Customized list of amplicons for cell selection
ampli.list <- read.csv(opt$ampli,header=FALSE)
sel.amplis <- as.character(ampli.list$V1)

# Apply cell selection function
data <- read.table(opt$file,
                     header = T)
thres <- resemble.cellfinder(data[,sel.amplis])
ful.cutoff <- data[,sel.amplis]>thres
ful.perc <- apply(ful.cutoff,1,function(x)sum(x)/length(x))
data <- data[ful.perc > val,]
# print number of cells and cut off used
print(paste("Cells detected:", nrow(data),", with cutoff:", val, ", using", length(sel.amplis), "amplicons"))
# save table with cells
write.table(data, opt$output,sep='\t', quote = FALSE)
