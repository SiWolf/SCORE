#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

create_gene_list <- function(sample_1){
  path_1 <- paste("../mapped/bowtie2/featureCounts/", sample_1, sep = "")
  setwd(path_1)
  gene_list = read.csv("counts", sep="", head=T, skip=1)[,c("Geneid")]
  return(gene_list)
}

create_count_matrix <- function(sample_2, sample_3, sample_4){
  counts_1 = read.csv("counts", sep="", head=T, skip=1, row.names = "Geneid")
  # Remove first six columns (genesymbol, chr, start, end, strand, length)
  counts_1 <- counts_1[ ,6:ncol(counts_1)]
  
  path_2 <- paste("../", sample_2, sep = "")
  setwd(path_2)
  counts_2 = read.csv("counts", sep="", head=T, skip=1, row.names = 1)
  counts_2 <- counts_2[ ,6:ncol(counts_2)]
  
  path_3 <- paste("../", sample_3, sep = "")
  setwd(path_3)
  counts_3 = read.csv("counts", sep="", head=T, skip=1, row.names = 1)
  counts_3 <- counts_3[ ,6:ncol(counts_3)]
  
  path_4 <- paste("../", sample_4, sep = "")
  setwd(path_4)
  counts_4 = read.csv("counts", sep="", head=T, skip=1, row.names = 1)
  counts_4 <- counts_4[ ,6:ncol(counts_4)]
  
  # Convert to matrix
  count_matrix <- cbind(counts_1, counts_2, counts_3, counts_4)
  count_matrix <- as.matrix(count_matrix)
  #head(count_matrix)
  return(count_matrix)
}

run_deseq2 <- function(read_counts, list_of_gene_names){
  library("DESeq2")
  
  # Assign condition
  condition <- factor(c(rep("normal", 2), rep("treated", 2)))
  
  # Create a coldata frame and instantiate the DESeqDataSet
  coldata <- data.frame(row.names=colnames(read_counts), condition)
  
  dds <- DESeqDataSetFromMatrix(countData=read_counts, colData=coldata, design=~condition)
  # Run the DESeq pipeline
  dds <- DESeq(dds)
  # Get differential expression results
  res <- results(dds)
  table(res$padj<0.05)
  # Order by adjusted p-value
  res <- res[order(res$padj), ]
  # Merge with normalized count data and gene symbols
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  new_resdata <- merge(as.data.frame(resdata), as.data.frame(list_of_gene_names), by="row.names", sort=FALSE)
  new_resdata <- new_resdata[3:13]
  names(new_resdata)[11] <- "Gene"
  head(resdata)
  
  # Write results
  path_final_1 = "../../../../"
  setwd(path_final_1)
  path_final_2 = "deg/"
  setwd(path_final_2)
  write.csv(new_resdata, file="deseq2-diffexpr-results.csv")
  # Plots
  hist(res$padj, breaks=50, col="grey")
  plotMA(res, ylim=c(-2,2))
  plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
}

run_edger <- function(Counts){
  library("edgeR")
  dgList <- DGEList(counts=Counts, genes=rownames(Counts))
  print("This shouldn't be running")
}

# Main

args<-commandArgs(TRUE)
argument_1 = args[1]
argument_2 = args[2]
argument_3 = args[3]
argument_4 = args[4]
argument_5 = args[5]

if (is.na(argument_1)){
  argument_1 = "DeSeq2"
  argument_2 = "SRR3199338"
  argument_3 = "SRR3199339"
  argument_4 = "SRR3199340"
  argument_5 = "SRR3199341"
}

gene_names <- create_gene_list(argument_2)
counts <- create_count_matrix(argument_3, argument_4, argument_5)

if (argument_1 == "DeSeq2"){
  run_deseq2(counts, gene_names)
} else{
  run_edger(counts)
}

#getwd()
#rld <- rlog(dds)