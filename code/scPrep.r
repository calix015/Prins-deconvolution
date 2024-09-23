####Header############################################################
# Project: Prins 2023: Identification of differential gene expression 
# responses in pulmonary cells from bulk RNA-Seq using deconvolution
# Code goal: scRNAseq in shape to run DWLS
# RI Analyst: Natalia Calixto Mancipe

## Libraries----------------------------------------------------------
.libPaths(c("/panfs/jay/groups/0/riss/calix015/myR_lib/",.libPaths()))
library(tidyverse)
library(future)
library(Seurat)
library(DWLS)
library(stringr)

## Setup--------------------------------------------------------------
# This session was run in OpenOnDemand using BigMem node 
# (32 core, 500Gb RAM) using R/4.1.0
setwd("~/2023/prins/")
# check the current active plan for parallelization with future
plan()
# change the current plan to access parallelization. 
# Parallelize: FindIntegrationAnchors, FindMarkers, and FindClusters 
plan("multicore", workers = 32)
# Need to increase "future.globals.maxSize" if using future R package 
# here I'm increasing to 200Gb
options(future.globals.maxSize = 200000 * 1024^2) # requires the number in bytes

outdir = "~/2023/prins/sc_data/"
if(!file.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Read ref data--------------------------------
# Read metadata
meta <- read_delim("sc_data/meta.tsv", delim = "\t",col_names = T, progress = T, 
                   show_col_types = T)
# This data is already transformed and given as ln(count+1)
# Just need to filter out the condition that is not of interest (SuHx)
holder <- read_delim("sc_data/exprMatrix.tsv", delim = "\t",col_names = T, 
                     progress = T, show_col_types = F)
# Get gene names
gnames <- holder$gene
holder <- data.matrix(holder[,2:ncol(holder)])
rownames(holder) <- gnames

# Create Seurat Object ---------
dataSC <- CreateSeuratObject(counts = holder, project = "scRNA_Hong") 
#head(dataSC@meta.data)
#str(dataSC)

# Add metadata (the first row is missing in the meta file provided)
nobs <- nrow(meta)
dataSC@meta.data <- dataSC@meta.data %>%
  mutate(barcode = meta$Cell[2:nobs],
         cell_type = meta$`cell type`[2:nobs],
         condition = meta$condition[2:nobs],
         ngene = meta$genes[2:nobs],
         nUMI = meta$UMIs[2:nobs],
         mt.pc = meta$`mito percent`[2:nobs])

Idents(dataSC) <- "condition"
#dim(dataSC) #20143 33391

# Get only the data we need -----
# Controls
scControl <- subset(dataSC, idents = "Control")
write.table(table(scControl$cell_type), 
            file = paste0(outdir,"scControl.cellnum.txt"),quote = F, sep = "\t",
            row.names = F, col.names = F)
saveRDS(scControl,file = paste0(outdir,"scControl.rds"))

# Monocrotaline
scMCT <- subset(dataSC, idents = "MCT")
write.table(table(scMCT$cell_type), 
            file = paste0(outdir,"scMCT.cellnum.txt"),quote = F, sep = "\t",
            row.names = F, col.names = F)
saveRDS(scControl,file = paste0(outdir,"scMCT.rds"))

# Visualize QC metrics as a violin plot ---------------------
# data is already filtered, so it should be good
pdf(paste0(outdir,"scControl.QC.pdf"), width = 7.5, height = 4.5)
VlnPlot(scControl, 
        features = c("ngene","nUMI","mt.pc","nCount_RNA"), 
        pt.size = 0, ncol = 4)
dev.off()

pdf(paste0(outdir,"scMCT.QC.pdf"), width = 7.5, height = 4.5)
VlnPlot(scMCT, 
        features = c("ngene","nUMI","mt.pc","nCount_RNA"), 
        pt.size = 0, ncol = 4)
dev.off()

## Cell type Markers in controls with DWLS - Seurat -----
# Remove spaces from cell names 
scControl@meta.data$c_type <- str_replace_all(scControl@meta.data$cell_type, pattern = "\\s+", replacement = "_")
head(scControl@meta.data)
Idents(scControl) <- "c_type"
c.type <- unique(Idents(scControl))
for(i in c.type){
  de_group <- FindMarkers(object=scControl, ident.1 = i, ident.2 = NULL, 
                          only.pos = TRUE, test.use = "bimod")
  saveRDS(de_group,file=paste0(outdir,"control/de_",i,".rds"))
}

## Cell type Markers in MCT with DWLS - Seurat -----
# Remove spaces from cell names 
scMCT@meta.data$c_type <- str_replace_all(scMCT@meta.data$cell_type, pattern = "\\s+", replacement = "_")
head(scMCT@meta.data)
Idents(scMCT) <- "c_type"
c.type <- unique(Idents(scMCT))
for(i in c.type){
  de_group <- FindMarkers(object=scMCT, ident.1 = i, ident.2 = NULL, 
                          only.pos = TRUE, test.use = "bimod")
  saveRDS(de_group,file=paste0(outdir,"MCT/de_",i,".rds"))
}

# Clean up
rm(list=c("outdir","de_group","i","counts"))

# SaveImage
save.image("~/2023/prins/scPrep.RData")