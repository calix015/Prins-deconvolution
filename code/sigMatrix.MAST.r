# Prins 2023
# Build Signature matrix using DEGs from single cell data

## Libraries----------------------------------------------------------
.libPaths(c("/panfs/jay/groups/0/riss/calix015/myR_lib/",.libPaths()))
library(tidyverse)
library(future)
library(Seurat)
library(DWLS)
library(Matrix)

## Setup--------------------------------------------------------------
# This session was run in OpenOnDemand using BigMem node 
# (32 core, 500Gb RAM) using R/4.1.0
setwd("~/2023/prins/sc_data/")
# Parallelize: FindIntegrationAnchors, FindMarkers, and FindClusters 
#plan("multicore", workers = 32)
#options(future.globals.maxSize = 200000 * 1024^2) # requires the number in bytes

load("~/2023/prins/scPrep.RData")

outdir = "~/2023/prins/sc_data/"
if(!file.exists(outdir)) dir.create(outdir, recursive = TRUE)

# SM of control cells -------------------------------------
#counts <- as.matrix(scControl@assays$RNA@counts)
#c.types <- Idents(scControl)
# Need to pass the full vector of cell ids (not just the unique ids!!)
#ModbuildSignatureMatrixMAST(scdata = counts, id = c.types, path = "/panfs/jay/groups/0/riss/calix015/2023/prins/sc_data/Control_MAST")

# SM of MCT cells -------------------------------------
#counts <- as.matrix(scMCT@assays$RNA@counts)
#c.types <- Idents(scMCT)
# Need to pass the full vector of cell ids (not just the unique ids!!)
#buildSignatureMatrixMAST(scdata = counts, id = c.types, path = "/panfs/jay/groups/0/riss/calix015/2023/prins/sc_data/MCT_MAST")


# Troubleshoot function -> Modified version! -----
# Gets errors, probably due to weird paths?
# Original function is 
ModbuildSignatureMatrixMAST<-function(scdata,id,path,diff.cutoff=0.5,
                                   pval.cutoff=0.01, f = 200)
  { #DEAnalysisMAST(scdata,id,path) this step is complete
  numberofGenes<-c()
  for (i in unique(id)){
    if(file.exists(paste(path,"/",i,"_MIST.RData", sep=""))){
      load(file=paste(path,"/",i,"_MIST.RData", sep=""))
      pvalue_adjusted<-p.adjust(cluster_lrTest.table[,3], method = "fdr",
                                n = length(cluster_lrTest.table[,3]))
      cluster_lrTest.table<-cbind(cluster_lrTest.table,pvalue_adjusted)
      DEGenes<-cluster_lrTest.table$Gene[intersect(which
                                              (pvalue_adjusted<pval.cutoff),
                      which(cluster_lrTest.table$log2fold_change>diff.cutoff))]
      nonMir = grep("MIR|Mir", DEGenes, invert = T)
      assign(paste("cluster_lrTest.table.",i,sep=""),
    cluster_lrTest.table[which(cluster_lrTest.table$Gene%in%DEGenes[nonMir]),])
      numberofGenes<-c(numberofGenes,length(DEGenes[nonMir]))
    }
  }

  #need to reduce number of genes
  #for each subset, order significant genes by decreasing fold change,
  #choose between 50 and 200 genes
  #for each, iterate and choose matrix with lowest condition number
  conditionNumbers<-c()
  for(G in 50:f){
    Genes<-c()
    j=1
    for (i in unique(id)){
      if(numberofGenes[j]>0){
        temp<-paste("cluster_lrTest.table.",i,sep="")
        print(temp)
        temp<-as.name(temp)
        temp<-eval(parse(text = temp))
        temp<-temp[order(temp$log2fold_change,decreasing=TRUE),]
        Genes<-c(Genes,varhandle::unfactor(temp$Gene[1:min(G,
                                                           numberofGenes[j])]))
      }
      j=j+1
    }
    Genes<-unique(Genes)
    #make signature matrix
    ExprSubset<-scdata[Genes,]
    Sig<-NULL
    for (i in unique(id)){
      Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
    }
    colnames(Sig)<-unique(id)
    conditionNumbers<-c(conditionNumbers,kappa(Sig))
  }
  G<-which.min(conditionNumbers)+min(49,numberofGenes-1)
  Genes<-c()
  j=1
  for (i in unique(id)){
    if(numberofGenes[j]>0){
      temp<-paste("cluster_lrTest.table.",i,sep="")
      print(paste0(temp,"2nd loop"))
      temp<-as.name(temp)
      temp<-eval(parse(text = temp))
      temp<-temp[order(temp$log2fold_change,decreasing=TRUE),]
      Genes<-c(Genes,varhandle::unfactor(temp$Gene[1:min(G,numberofGenes[j])]))
    }
    j=j+1
  }
  Genes<-unique(Genes)
  ExprSubset<-scdata[Genes,]
  Sig<-NULL
  for (i in unique(id)){
    Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig)<-unique(id)
  save(Sig,file=paste(path,"/Sig.RData",sep=""))
  saveRDS(Sig,file=paste(path,"/Sig.rds",sep=""))
  return(Sig)
}

# Problem with "-" symbol in a name
# Fix for control dataset
scControl@meta.data <- scControl@meta.data %>%
  mutate(c_type = case_when(c_type == "Non-classical_monocytes" ~ "Non_classical_monocytes", T ~ c_type))
Idents(scControl) <- "c_type"
counts <- as.matrix(scControl@assays$RNA@counts)
c.types <- Idents(scControl)
path <- "/panfs/jay/groups/0/riss/calix015/2023/prins/sc_data/Control_MAST"
ModbuildSignatureMatrixMAST(scdata = counts, id = c.types, path = "/panfs/jay/groups/0/riss/calix015/2023/prins/sc_data/Control_MAST")

# Fix for MCT dataset
scMCT@meta.data <- scMCT@meta.data %>%
  mutate(c_type = case_when(c_type == "Non-classical_monocytes" ~ "Non_classical_monocytes", T ~ c_type))
Idents(scMCT) <- "c_type"
counts <- as.matrix(scMCT@assays$RNA@counts)
c.types <- Idents(scMCT)
path <- "/panfs/jay/groups/0/riss/calix015/2023/prins/sc_data/MCT_MAST/"
ModbuildSignatureMatrixMAST(scdata = counts, id = c.types, path = "/panfs/jay/groups/0/riss/calix015/2023/prins/sc_data/MCT_MAST")

save.image("~/2023/prins/scPrep.RData")