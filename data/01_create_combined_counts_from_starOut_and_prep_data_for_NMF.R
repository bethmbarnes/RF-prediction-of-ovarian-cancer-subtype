#load libraries
library(stringr)
library(DESeq2)
library(dplyr)

#set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/'))
setwd(paste0(myDir[-length(myDir)], collapse = '/'))
getwd()

#Import read counts
files <- dir(path="./data/starOut", full.names=T, pattern="ReadsPerGene.out.tab")
counts <- c()
for( i in seq_along(files) ){
  x <- read.table(file=files[i], sep="\t", header=F, as.is=T, skip=4)
  counts <- cbind(counts, x[,2])
}

#tidy row and column names
files_colnames <- str_replace(files, pattern = "ReadsPerGene.out.tab", replacement = "")
files_colnames <- str_replace(files_colnames, pattern = "./data/starOut/", replacement = "")
colnames(counts) <- files_colnames
rownames(counts) <- x$V1

#write to CSV
#write.csv(counts, file="data/combined_counts_realigned.csv")
counts <- read.csv(file="data/combined_counts_realigned.csv", row.names=1)

#targets file
targets <- read.csv(file="data/CCLE_targets.csv", row.names=1)

#check order and naming
all(rownames(targets) %in% colnames(counts)) 
all(rownames(targets) == colnames(counts))

#make cell line name the row/column name
colnames(counts) <- targets$Cell_Line
rownames(targets) <- targets$Cell_Line

#remove ENSG version number
rownames(counts) <- str_replace(rownames(counts), pattern = ".[0-9]+$", replacement = "")

#get gene symbols from biomart (output provided in data)
#ENSG IDs to use
#genes <-  rownames(counts)
#mart to use
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#run biomart query. Takes ~5 min
#G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"), values = genes, mart = mart)
#write.csv(G_list, file="data/biomart_ENSG_to_HGNC_symbol.csv")
G_list <- read.csv(file="data/biomart_ENSG_to_HGNC_symbol.csv")

#merge counts with biomart return
counts_hgnc <- merge(G_list, counts, by.x="ensembl_gene_id", by.y=0)
counts_hgnc <- counts_hgnc[,-1:-2]
counts_hgnc <- as.data.frame(counts_hgnc %>%
                               group_by(hgnc_symbol) %>% 
                               summarise_all(funs(sum)))
rownames(counts_hgnc) <-counts_hgnc[,1]
counts_hgnc <- counts_hgnc[-1,-1]

#VST transformation using DESeq2 (output provided in data)
dds <- DESeqDataSetFromMatrix(countData = counts_hgnc,
                              colData = targets,
                              design = ~ 1) #1 makes it blind to column data
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vsd <- as.data.frame(vsd@assays@data@listData[[1]])
#write.csv(vsd, file="data/combined_counts_realigned_variance_stabilised_transformed.csv")