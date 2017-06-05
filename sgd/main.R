rm(list=ls())

#Parameter 
medium <- "HAM" #HAM or FBS

##In silico
insilicoSGD.files <- dir("predictions/")[c(grep("TCGA",dir("predictions/")),grep("HMR",dir("predictions/")),grep("Recon3",dir("predictions/")))]
insilicoSGD.files <- insilicoSGD.files[grepl(medium,insilicoSGD.files)]

#Get all genes tested
converttoHGNC <- function(genes,from){
  genes <- unlist(lapply(genes,function(x) strsplit(x,"\\.")[[1]][1]))
  if (from=="ensembl"){
    load("~/Documents/Data/Geneannot/annotation_ensembl_gene_id.RData")
    genes <- annotEntrez[genes,"hgnc_symbol"]
  } else if (from=="entrez"){
    load("~/Documents/Data/Geneannot/annotation_entrezgene.RData")
    genes <- annotEntrez[genes,"hgnc_symbol"]
  } else {stop("from must be either 'ensembl' or 'entrez'")}
  genes
}

genes   <- c()
for (file in insilicoSGD.files){
  sgd   <- readLines(paste0("predictions/",file))
  if (all(grepl("ENSG",sgd))){
    sgd <- converttoHGNC(sgd,"ensembl")
  } else if (all(grepl("\\.",sgd))){
    sgd <- converttoHGNC(sgd,"entrez")
  } else {
    stop("All genes must be either ensembl or entrez in ",file)
  }
  genes <- union(genes,sgd)
}


#Compile matrix of SGD
sgdmat   <- matrix(0,nrow=length(genes),ncol=length(insilicoSGD.files),
                   dimnames=list(genes,insilicoSGD.files))
for (file in insilicoSGD.files){
  sgd   <- readLines(paste0("predictions/",file))
  if (all(grepl("ENSG",sgd))){
    sgd <- converttoHGNC(sgd,"ensembl")
  } else if (all(grepl("\\.",sgd))){
    sgd <- converttoHGNC(sgd,"entrez")
  } else {
    stop("All genes must be either ensembl or entrez in ",file)
  }
  sgdmat[sgd,file] <- 1
}

#Reannotate file names
ID      <- c()
for (file in insilicoSGD.files){
  if (grepl("HMR",file)&grepl("FBS",file)) {
    ID <- c(ID,"HMR2_FBS")
  } else if (grepl("HMR",file)&grepl("HAM",file)) {
    ID <- c(ID,"HMR2_HAM")
  } else if (grepl("Recon3",file)&grepl("FBS",file)) {
    ID <- c(ID,"Recon3_FBS")
  } else if (grepl("Recon3",file)&grepl("HAM",file)) {
    ID <- c(ID,"Recon3_HAM")
  } else {
    ID <- c(ID,paste0(substr(file,14,29),"_HAM"))
  }
}
colnames(sgdmat) <- ID

#Heatmap
library(gplots)
library(RColorBrewer)
save    <- T
matrix  <- sgdmat
hmcols  <- colorRampPalette(brewer.pal(9,"Greys"))(2)
label   <- rep("GBM",ncol(sgdmat))
  label[grep("HMR",colnames(sgdmat))] <- "HMR2"
  label[grep("Recon",colnames(sgdmat))] <- "Recon3"
labcols <- brewer.pal(nlevels(factor(label)),"Set2")[factor(label)]
if (save) pdf(file = paste0("plots/SGD_heatmap_GBMvsHMR2vsRecon3_",medium,".pdf"), width=4, height=8)
hm <- heatmap.2(matrix,scale="none",col=hmcols,
                density.info="none",cexCol = 0.8,symm=F,
                labRow=toupper(rownames(matrix)),labCol=toupper(colnames(matrix)),
                trace="none",ColSideColors=labcols)
legend("topright",legend = toupper(levels(factor(label))),
       col=brewer.pal(nlevels(factor(label)),"Set2"), 
       lty= 1, lwd = 5,cex=0.8)
if (save) dev.off()
