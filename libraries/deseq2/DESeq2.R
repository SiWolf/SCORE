#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

library("DESeq2")

args<-commandArgs(TRUE)
a1 = args[1]
a2 = args[2]
a3 = args[3]
a4 = args[4]

#a1 = "SRR3199338"
#a2 = "SRR3199339"
#a3 = "SRR3199340"
#a4 = "SRR3199341"

#getwd()

path_1 <- paste("../../mapped/bowtie2/featureCounts/", a1, sep = "")
setwd(path_1)
gene_names = read.csv("counts", sep="", head=T, skip=1)[,c("Geneid")]
counts_1 = read.csv("counts", sep="", head=T, skip=1, row.names = "Geneid")
# Remove first six columns (genesymbol, chr, start, end, strand, length)
counts_1 <- counts_1[ ,6:ncol(counts_1)]

path_2 <- paste("../", a2, sep = "")
setwd(path_2)
counts_2 = read.csv("counts", sep="", head=T, skip=1, row.names = 1)
counts_2 <- counts_2[ ,6:ncol(counts_2)]

path_3 <- paste("../", a3, sep = "")
setwd(path_3)
counts_3 = read.csv("counts", sep="", head=T, skip=1, row.names = 1)
counts_3 <- counts_3[ ,6:ncol(counts_3)]

path_4 <- paste("../", a4, sep = "")
setwd(path_4)
counts_4 = read.csv("counts", sep="", head=T, skip=1, row.names = 1)
counts_4 <- counts_4[ ,6:ncol(counts_4)]

# Convert to matrix
counts <- cbind(counts_1,counts_2,counts_3,counts_4)
counts <- as.matrix(counts)
head(counts)

# Assign condition
condition <- factor(c(rep("normal", 2), rep("treated", 2)))

# Create a coldata frame and instantiate the DESeqDataSet
coldata <- data.frame(row.names=colnames(counts), condition)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
# Run the DESeq pipeline
dds <- DESeq(dds)
# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
# Order by adjusted p-value
res <- res[order(res$padj), ]
# Merge with normalized count data and gene symbols
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
new_resdata <- merge(as.data.frame(resdata), as.data.frame(gene_names), by="row.names", sort=FALSE)
new_resdata <- new_resdata[3:13]
names(new_resdata)[11] <- "Gene"
head(resdata)
# Write results
path_final_1 = "../../../../"
setwd(path_final_1)
path_final_2 = "deg/"
setwd(path_final_2)
write.csv(new_resdata, file="deseq2-diffexpr-results.csv")
# Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")
# After multiple testing correction
hist(res$padj, breaks=50, col="grey")

# Removing N/As?
#results(dds, independentFiltering=FALSE)
#res$pvalue[res$baseMean < x] <- NA
#res$padj <- p.adjust(res$pvalue, method="BH")