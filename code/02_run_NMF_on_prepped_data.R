#load libraries
library(NMF)
library(Biobase)

#set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/'))
setwd(paste0(myDir[-length(myDir)], collapse = '/'))
getwd()

#targets file
targets <- read.csv(file="data/CCLE_targets.csv", row.names=2)

#Variance stabilsied reads file
vsd <- read.csv(file="data/combined_counts_realigned_variance_stabilised_transformed.csv", row.names=1)

#Calculate median absolute deviation
vsd$rmads <- apply(vsd, 1, mad)
#subset vsd and remove MAD column
vsd_sub <-subset(vsd[order(vsd$rmads),], rmads > 1.5, select=1:44)

#column data meta info for eset
metadata <- data.frame(labelDescription= c("Annotated subtype", "Domcke score"),
                       row.names=c("Annotated_Subtype", "Domcke_Score"))

targets_temp <- as.data.frame(targets[,-1:-2])
rownames(targets_temp) <- rownames(targets)

phenoData <- new("AnnotatedDataFrame", data=targets_temp, varMetadata=metadata)

eSet <- ExpressionSet(assayData=(as.matrix(vsd_sub)), phenoData=phenoData)

#NMF rank test for k value of 2 to 10
estim.r <- nmf(eSet, 2:10, nrun=5, seed=123456) #nrun=50
#view rank survey
plot(estim.r)
#view all consensus maps
consensusmap(estim.r, annCol=eSet, labCol=NA, labRow=NA)

#NMF for k of 5
estim.r.5 <- nmf(eSet, 5, nrun=20) #nrun=200
#view consensus map
consensusmap(estim.r.5, labRow=NA, annCol=eSet)
#view basis map
basismap(estim.r.5, subsetRow=TRUE, labRow=NA)
#extract silhuette scores (gives assignment of cell lines)
si <- silhouette(estim.r.5, what = 'consensus')
#write.csv(si, file="silhuette_scores.csv")
#Extract metagenes
s <- extractFeatures(estim.r.5) 
dat <- eSet@assayData[["exprs"]]
#get indices of metagenes
basis_k5_1<- s[[1]]
basis_k5_2<- s[[2]]
basis_k5_3<- s[[3]]
basis_k5_4<- s[[4]]
basis_k5_5<- s[[5]]
#use indices to extract data
genes_k5_1<-dat[basis_k5_1,]
genes_k5_2<-dat[basis_k5_2,]
genes_k5_3<-dat[basis_k5_3,]
genes_k5_4<-dat[basis_k5_4,]
genes_k5_5<-dat[basis_k5_5,]
genes_k5_all <- rbind(genes_k5_1, genes_k5_2, genes_k5_3, genes_k5_4, genes_k5_5)
#write.csv(genes_k5_all, file="data/NMF_metagenes_k5_all.csv")

#NMF for k of 2 (from supplementary figures)
estim.r.2 <- nmf(eSet, 2, nrun=20) #nrun=200
#view consensus map
consensusmap(estim.r.2, labRow=NA, annCol=eSet)
#view basis map
basismap(estim.r.2, subsetRow=TRUE, labRow=NA)
#extract silhuette scores (gives assignment of cell lines)
si <- silhouette(estim.r.2, what = 'consensus')
#write.csv(si, file="silhuette_scores.csv")
#Extract metagenes
s <- extractFeatures(estim.r.2) 
dat <- eSet@assayData[["exprs"]]
#get indices of metagenes
basis_k2_1<- s[[1]]
basis_k2_2<- s[[2]]
#use indices to extract data
genes_k2_1<-dat[basis_k2_1,]
genes_k2_2<-dat[basis_k2_2,]
genes_k2_all <- rbind(genes_k2_1, genes_k2_2)
#write.csv(genes_k2_all, file="NMF_metagenes_k2_all.csv")