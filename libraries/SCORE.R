#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

library("DESeq2")

create_gene_list <- function(sample){
  sample_path <- paste("mapped/bowtie2/featureCounts/", sample, sep = "")
  setwd(sample_path)
  gene_list = read.csv("counts", sep="", head=T, skip=1)[,c("Geneid")]
  return(gene_list)
}

create_count_matrix <- function(sample_list, gene_list){
  sample_nr = 0
  # Lets receive the rest of the counts...
  for (sample in sample_list){
    sample_nr = sample_nr + 1
    path <- paste("../", sample, sep = "")
    setwd(path)
    counts = read.csv("counts", sep="", head=T, skip=1, row.names = 1)
    # Remove first six columns (genesymbol, chr, start, end, strand, length)
    counts <- counts[ ,6:ncol(counts)]
    if (sample_nr == 1){
      count_matrix <- counts
    } else {
      count_matrix <- cbind(count_matrix, counts)  
    }
  }
  
  # Convert to matrix
  count_matrix <- as.matrix(count_matrix)
  # Rename columns
  colnames(count_matrix) <- paste(sample_list)
  rownames(count_matrix) <- paste(gene_list)
  # head(count_matrix)
  return(count_matrix)
}

run_deseq2 <- function(list_of_gene_names, sample_counts, sample_conditions){
  # Assign condition
  # condition <- factor(c(rep("normal", 2), rep("treated", 2)))
  
  # Create a coldata frame and instantiate the DESeqDataSet
  coldata <- data.frame(row.names=colnames(sample_counts), sample_conditions)
  
  dds <- DESeqDataSetFromMatrix(countData=sample_counts, colData=coldata, design=~sample_conditions)
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
  results_folder = "../../../../deg/"
  setwd(results_folder)
  write.csv(new_resdata, file="deseq2-diffexpr-results.csv")
  # Plots
  hist(res$padj, breaks=50, col="grey")
  plotMA(res, ylim=c(-2,2))
  plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
}

run_edger <- function(read_counts){
  library("edgeR")
  dgList <- DGEList(counts=read_counts, genes=rownames(read_counts))
  print("This shouldn't be running")
}

# Main
args<-commandArgs(TRUE)
argument_1 = args[1]
argument_2 = args[2]

if (is.na(argument_1)){
  argument_1 = "DeSeq2"
  argument_2 = "Metadata.tsv"
  setwd("../")
}

metadata = read.table(file = paste("raw/", argument_2, sep = ""), sep = "\t", header = FALSE)
# metadata = read.table(file="raw/Metadata.tsv", sep = "\t", header = FALSE)
gene_names <- create_gene_list(metadata$V1[1])
gene_counts <- create_count_matrix(metadata$V1, gene_names)

if (argument_1 == "DeSeq2"){
  run_deseq2(gene_names, gene_counts, metadata$V2)
} else{
  run_edger(gene_counts)
}
