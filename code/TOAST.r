####Header############################################################
# Project: Prins 2023: Identification of differential gene expression 
# responses in pulmonary cells from bulk RNA-Seq using deconvolution
# Code goal: DEG on deconvoluted samples
# RI Analyst: Natalia Calixto Mancipe

## Libraries----------------------------------------------------------
.libPaths(c("/panfs/jay/groups/0/riss/calix015/myR_lib/",.libPaths()))
library(tidyverse)
library("matrixStats")
#library(pheatmap)
library(TOAST)
#library(org.Rn.eg.db)
library(clusterProfiler)
#library(edgeR)

## Setup--------------------------------------------------------------
# This session was run in OpenOnDemand using 2 cores and 1 node, R/4.1.0
# Vignette analysis code
# https://raw.githubusercontent.com/cozygene/TCA/master/vignettes/vignette_analysis.R

setwd("~/2023/prins/")
load("./toast.Rdata")
#indir="~/2023/prins/bulk_data/"
outdir = "~/2023/prins/degs/GSEA/"
if(!file.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Functions ------------------------------------------------------------------
gsea_plots <- function(gsea_holder, cell_type, test_type){
  
  isets <- rownames(gsea_holder@result %>%
                    select(p.adjust, NES) %>%
                    filter(NES<0 & p.adjust < 0.05))
  for (i in isets){
    plot2save <- gseaplot(gsea_holder, geneSetID = i, by = "runningScore",title = i)
    plotfile <- paste0(outdir,cell_type,test_type,i,".png")
    ggsave(filename = plotfile, plot = plot2save, width = 5, height = 3)
  }
}

## Prepare GSEA Universe (SKIP) ----------------------------------------------
library(msigdbr)
msigdbr_species() # use Rattus norvegicus
m_t2g <- msigdbr(species = "Rattus norvegicus", category = "H") %>% 
  dplyr::select(gs_name, ensembl_gene)
head(m_t2g)
dim(m_t2g) #7300    2
# Filter out genes that were not detected for universe
m_t2g <- m_t2g %>%
  filter(ensembl_gene %in% rownames(reads))
dim(m_t2g) #6683    2

## Prepare bulk data SKIP -----
reads <- read_delim(paste0(indir,"subread_counts.txt"), delim = "\t", comment = "#")
colnames(reads)

# Should be nicely put for differential exression
# Only using genes > 300 bp
reads <- reads %>%
  filter(Length > 300) %>%
  select(-all_of(c("Chr","Start","End","Length","Strand")))


meta <- tibble(sample = colnames(reads)[2:dim(reads)[2]],
               treatment = c(rep("control",4),rep("MCT_ferrostatin",4),rep("MCT",4)),
               ferrostatin = c(rep("control",4),rep("ferrostatin",4),rep("vehicle",4)),
               MCT = c(rep("control",4),rep("mct",8)))

genes <- reads[,1]
reads <- data.matrix(reads[,2:ncol(reads)])
rownames(reads) <- genes$Geneid
dim(reads) #25680    12

# DGEList, filter by expression
# grouping by treatment and get rid of lowly expressed genes
y <- DGEList(counts = reads, group = meta$treatment, genes = genes$Geneid)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

# Get the normalized counts using cpm
# nc <- cpm(y)
# rownames(nc) <- y$genes$Geneid #nc contains the normalized reads
reads <- y$counts
dim(reads) #15693    12
rm(list=c("y","keep"))

## Get composition matrix (SKIP) ----------------------------------------------
# Proportion data and make df, then numeric matrix
file_names <- list.files(outdir, pattern = "^DWLS_NV_", full.names = T)
Prop <- NULL
for (i in 1:length(file_names)){
  x.file <- file_names[i]
  test <- str_split(x.file, pattern = "/", simplify = T)
  test <- test[,dim(test)[2]]
  test <- str_replace(test, pattern = ".txt",replacement = "")
  test <- str_replace(test, pattern = "DWLS_",replacement = "")
  
  a <- read_delim(x.file, delim = "\t", col_names = F) 
  colnames(a) <- c("cell_type",test)
  if(i == 1){
    Prop <- a
  }
  if(i != 1){
    Prop <- cbind(Prop, a[,2])
  }
}
rm(list = c("a","file_names","i","test","x.file"))

#Heatmaps of devonvolution data 
cell_cluster <- Prop$cell_type
ct.percent <- Prop %>%
  pivot_longer(-cell_type,names_to = "sample", values_to = "prop") %>%
  mutate(ct.pc = prop*100) %>%
  pivot_wider(id_cols = -prop, names_from = sample, values_from = ct.pc) %>%
  select(-cell_type)
ct.percent <- data.matrix(ct.percent)
rownames(ct.percent) <- cell_cluster

annot_col <- data.frame(treatment=meta$treatment)
rownames(annot_col)<- meta$sample

# Sample composition
pdf("./deconv/heatmap.pc.pdf")
pheatmap(ct.percent, main = "Sample cell-type composition (%)", annotation_col = annot_col, show_colnames = F, annotation_names_col = F)
dev.off()

# Correlation of cell types based on predicted abundance
pdf("deconv/ct.cor.pdf")
pheatmap(cor(t(ct.percent)), main = "Pearson Correlation of Cell Types")
dev.off()

# Correlation of samples based on cell type abundance
p1 <- pheatmap(cor(ct.percent), main = "Pearson Correlation of Samples",annotation_col = annot_col, show_colnames = F)
pdf("deconv/sample.cor.pdf")
p1
dev.off()

## Group cell types -----
# Interesting cell types for researcher: Endothelial, Macrophages Smooth muscle
ct.prop <- Prop %>%
  pivot_longer(-cell_type,names_to = "sample", values_to = "prop") %>%
  mutate(cell_cluster = case_when(cell_type %in% c("Alveolar_macrophages","Proliferating_macrophages","Interstitial_macrophages") ~ "Macrophages",
                            cell_type %in% c("AT2_cells","AT1_cells") ~ "other",
                            cell_type %in% c("Endothelial_arterial_1","Endothelial_arterial_2","Endothelial_capillary") ~ "Endothelial",
                            cell_type %in% c("Neutrophils","NK_cells_1","NK_cells_2","Naive_T_cells",
                                             "Regulatory_T_cells", "Proliferating_T_cells","CD8_T_cells",
                                             "T_cells_Serpinb6", "B_cells","ILC2","Ciliated_cells",
                                             "Non_classical_monocytes","Classical_monocytes","Conventional_dendritic",
                                             "Plasmacytoid_dendritic", "Club_cells", "Fibroblasts","Mesothelial_cells") ~ "other",
                            cell_type %in% c("Smooth_muscle_cells") ~ "other",
                            cell_type %in% c("Mast_cells") ~ "other",
                            T ~ "other")) %>%
  ungroup()%>%
  group_by(pick(sample,cell_cluster)) %>%
  summarise(comb_prop = sum(prop)) %>%
  pivot_wider(names_from = sample, values_from = comb_prop)

# Cant deconvolve the control vs MCT samples using the grouping of cells 
# (Endothelial, Macrophages, smooth muscle, other) results are NA

## Prepare composition matrix Endothelial, Macrophage and other (SKIP)---------
cell_cluster <- ct.prop$cell_cluster
ct.prop <- ct.prop%>%select(-cell_cluster)
#reorder coulumns
ct.prop <- ct.prop %>% select(all_of(colnames(reads)))
# Convert to matrix
ct.prop <- data.matrix(ct.prop)
# add cell names
rownames(ct.prop) <- cell_cluster
# transpose 
ct.prop <- t(ct.prop) # ct.prop has the cell composition

# Check sample correlation with chosen cell grouping
p1 <- pheatmap(cor(ct.prop), main = "Cell type group abundance correlation")
png("./degs/Macr_end_other.cor.png")
p1
dev.off()

p2 <- pheatmap(cor(t(ct.prop)), main = "Sample correlation (Macrophages, Endothelial, other)", show_colnames = F)
png("./degs/Macr_end_other.sample_cor.png")
p2
dev.off()

## Test #1: Control vs MCT ----------------------------------------------------
#subset matrices
samples <- meta %>% filter(treatment != "MCT_ferrostatin")
reads.subset <- reads[,samples$sample] # use raw reads
head(reads.subset)
dim(reads.subset) #15693     8
ct.prop.subset <- ct.prop[samples$sample,]

# Subset metadata
meta.subset <- meta %>%
  filter(treatment != "MCT_ferrostatin")
annot_col <- data.frame(treatment = factor(meta.subset$treatment))
rownames(annot_col)<- meta.subset$sample

# Fit model
Design_out <- makeDesign(annot_col, ct.prop.subset)
fitted_model <- fitModel(Design_out, reads.subset)
fitted_model$all_coefs # list all phenotype names
fitted_model$all_cell_types # list all cell type names

# Test Macrophages
macrph_table <- csTest(fitted_model, coef = "treatment",cell_type = "Macrophages")
head(macrph_table)
summary(macrph_table$fdr)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9998  0.9998  0.9998  0.9998  0.9998  1.0000 
# macrph_table %>%
#   arrange(fdr) %>%
#   write.table(file="Macrophages_control_MCT_degs.txt",quote = F, sep = "\t",row.names = T, col.names = T)

# Test endothelial
end_table <- csTest(fitted_model, coef = "treatment", cell_type = "Endothelial")
head(end_table)
summary(end_table$fdr)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1       1       1       1       1       1 
# end_table %>%
#   arrange(fdr) %>%
#   write.table(file="Endothelial_control_MCT_degs.txt",quote = F, sep = "\t",row.names = T, col.names = T)

# Test the effect of treatment in all cell types
res_table <- csTest(fitted_model, coef = "treatment", cell_type = "joint")
head(res_table)
summary(res_table$fdr)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3846  0.4155  0.5110  0.5697  0.6955  0.9993 
# res_table %>%
#   arrange(fdr) %>%
#   write.table(file="Joint_control_MCT_degs.txt",quote = F, sep = "\t",row.names = T, col.names = T)

## GSEA analysis of DEGs in Control vs MCT --------------------------------------
# Organized by decreasing FDR, so only the genes in the bottom of the list 
# are sort of differentially expressed between conditions
# So only the significant results with negative NES are based on DE

# Macrophages
gene_list <- macrph_table$fdr
names(gene_list) <- rownames(macrph_table)
gene_list <- sort(gene_list, decreasing = T)
mac_gsea <- GSEA(geneList = gene_list, TERM2GENE = m_t2g, exponent = 1, nPermSimple = 100000)  
saveRDS(mac_gsea, paste0(outdir,"MCT_control.Macrophages.gsea.rds"))
mac_gsea@result %>% write.table(file="./degs/GSEA/Macrophages_control_MCT_GSEA.txt",quote = F, sep = "\t",row.names = T, col.names = T)
gsea_plots(mac_gsea, "Macrophages_", "control_MCT_")

# Endothelial
gene_list <- end_table$fdr
names(gene_list) <- rownames(end_table)
gene_list <- sort(gene_list, decreasing = T)
end_gsea <- GSEA(geneList = gene_list, TERM2GENE = m_t2g, exponent = 1, nPermSimple = 100000)  
saveRDS(end_gsea, paste0(outdir,"MCT_control.Endothelial.gsea.rds"))
end_gsea@result %>% write.table(file="./degs/GSEA/Endothelial_control_MCT_GSEA.txt",quote = F, sep = "\t",row.names = T, col.names = T)
gsea_plots(end_gsea, "Endothelial_", "control_MCT_")

# Joint
gene_list <- res_table$fdr
names(gene_list) <- rownames(res_table)
gene_list <- sort(gene_list, decreasing = T)
joint_gsea <- GSEA(geneList = gene_list, TERM2GENE = m_t2g, exponent = 1, nPermSimple = 100000)  
saveRDS(joint_gsea, paste0(outdir,"MCT_control.all_cells.gsea.rds"))
joint_gsea@result %>% write.table(file="./degs/GSEA/Joint_control_MCT_GSEA.txt",
                                  quote = F, sep = "\t",row.names = T, col.names = T)
gsea_plots(joint_gsea, "all_cells_", "control_MCT_")

## Test #2: MCT vs MCT_Ferrostatin-------------------------------------------------
#subset matrices
samples <- meta %>% filter(treatment != "control")
reads.subset <- reads[,samples$sample] 
head(reads.subset)
dim(reads.subset) #15693     8
ct.prop.subset <- ct.prop[samples$sample,]

# Subset metadata
meta.subset <- meta %>%
  filter(treatment != "control")
annot_col <- data.frame(treatment = factor(meta.subset$treatment))
rownames(annot_col)<- meta.subset$sample

# Fit model
Design_out <- makeDesign(annot_col, ct.prop.subset)
fitted_model <- fitModel(Design_out, reads.subset)
fitted_model$all_coefs 
fitted_model$all_cell_types

# Test Macrophages
macrph_table <- csTest(fitted_model, coef = "treatment",cell_type = "Macrophages")
head(macrph_table)
summary(macrph_table$fdr)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.7630  0.7630  0.7630  0.7697  0.7630  0.9999 
# macrph_table %>%
#   arrange(fdr) %>%
#   write.table(file="Macrophages_MCT_fer_degs.txt",quote = F, sep = "\t",row.names = T, col.names = T)

# Test endothelial
end_table <- csTest(fitted_model, coef = "treatment", cell_type = "Endothelial")
head(end_table)
summary(end_table$fdr)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.8668  0.8668  0.8668  0.8721  0.8668  1.0000 
# end_table %>%
#   arrange(fdr) %>%
#   write.table(file="Endothelial_MCT_fer_degs.txt",quote = F, sep = "\t",row.names = T, col.names = T)

# Test the effect of treatment in all cell types
res_table <- csTest(fitted_model, coef = "treatment", cell_type = "joint")
head(res_table)
summary(res_table$fdr)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1       1       1       1       1       1 
# res_table %>%
#   arrange(fdr) %>%
#   write.table(file="Joint_MCT_fer_degs.txt",quote = F, sep = "\t",row.names = T, col.names = T)

## GSEA analysis of DEGs in MCT vs MCT_Ferrostatin ----------------------------
# Macrophages
gene_list <- macrph_table$fdr
names(gene_list) <- rownames(macrph_table)
gene_list <- sort(gene_list, decreasing = T)
mac_gsea <- GSEA(geneList = gene_list, TERM2GENE = m_t2g, exponent = 1, nPermSimple = 100000)  
saveRDS(mac_gsea, paste0(outdir,"MCT_drug.Macrophages.gsea.rds"))
mac_gsea@result %>% write.table(file="./degs/GSEA/Macrophages_MCT_fer_GSEA.txt",quote = F, sep = "\t",row.names = T, col.names = T)
gsea_plots(mac_gsea, "ferrostatin_MCT/Macrophages_", "drug_MCT_")

# Endothelial
gene_list <- end_table$fdr
names(gene_list) <- rownames(end_table)
gene_list <- sort(gene_list, decreasing = T)
end_gsea <- GSEA(geneList = gene_list, TERM2GENE = m_t2g, exponent = 1, nPermSimple = 100000)  
saveRDS(end_gsea, paste0(outdir,"MCT_drug.Endothelial.gsea.rds"))
end_gsea@result %>% write.table(file="./degs/GSEA/Endothelial_MCT_fer_GSEA.txt",quote = F, sep = "\t",row.names = T, col.names = T)
gsea_plots(end_gsea, "ferrostatin_MCT/Endothelial_", "drug_MCT_")

# Joint
gene_list <- res_table$fdr
names(gene_list) <- rownames(res_table)
gene_list <- sort(gene_list, decreasing = T)
joint_gsea <- GSEA(geneList = gene_list, TERM2GENE = m_t2g, exponent = 1, nPermSimple = 100000)  
saveRDS(joint_gsea, paste0(outdir,"MCT_drug.all_cells.gsea.rds"))
joint_gsea@result %>% write.table(file="Joint_MCT_fer_GSEA.txt",
                                  quote = F, sep = "\t",row.names = T, col.names = T)
gsea_plots(joint_gsea, "ferrostatin_MCT/all_cells_", "drug_MCT_")

## GSEA figures -----
setwd("~/2023/prins/degs/GSEA/")
outdir = "~/2023/prins/degs/GSEA/"
if(!file.exists(outdir)) dir.create(outdir, recursive = TRUE)
library(RColorBrewer)

description <- NULL
NES <- NULL
FDR <- NULL
cells <- NULL
tests <- NULL
gene_ratios <- NULL

file_names <- list.files(pattern = "_GSEA.txt", full.names = T)
for(i in file_names){
  file.name <- str_split(i, pattern = "/", simplify = T)
  file.name <- file.name[length(file.name)]
  file.name <- str_remove(file.name, pattern = "_GSEA.txt")
  cell <- str_split(file.name, pattern = "_", simplify = T)[1]
  test <- paste(str_split(file.name, pattern = "_", simplify = T)[2],
                str_split(file.name, pattern = "_", simplify = T)[3],
                sep = ".vs.")
  
  gsea <- read_delim(i)
  gene.ratio <- str_count(gsea$core_enrichment, pattern = "/")/gsea$setSize
  
  gene_ratios <- append(gene_ratios, gene.ratio)
  description <- append(description, gsea$Description)
  NES <- append(NES, gsea$NES)
  FDR <- append(FDR, gsea$p.adjust)
  cells <- append(cells, rep(cell, length(gsea$NES)))
  tests <- append(tests, rep(test,length(gsea$NES)))
}

gsea2plot <- tibble(cells, tests, description, NES, FDR, gene_ratios)
colnames(gsea2plot) <- c("cells","tests","description","NES","FDR","Gene_ratio")

gsea2plot <- gsea2plot %>%
  filter(NES < 0 & FDR < 0.05) %>%
  mutate(NES = -1*NES,
         cells = factor(cells, levels = c("Endothelial","Macrophages","Joint"),
                        labels = c("Endothelial","Macrophages","All_cells")),
         pathway= str_remove(description, pattern = "HALLMARK_"),
         pathway= str_replace_all(pathway, pattern = "_", replacement = " "))

gsea2plot %>%
  write.table(file="GSEA.summary.txt",quote = F, sep = "\t",row.names = F, col.names = T)

unwanted <- c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
              "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
              "HALLMARK_G2M_CHECKPOINT",
              "HALLMARK_E2F_TARGETS",
              "HALLMARK_DNA_REPAIR")

holder <- gsea2plot %>% 
  filter(tests == "MCT.vs.fer") %>%
  mutate(keep = case_when(description %in% unwanted ~ "no", 
                          T ~ "yes")) %>%
  filter(keep == "yes")

ggplot(holder, aes(x=cells, y=pathway))+
  geom_point(aes(size = Gene_ratio, col=FDR))+
  ylab(element_blank())+
  xlab(element_blank())+
  scale_color_viridis_c()+
  ggtitle("Pathways affected by ferrostatin in MCT individuals", subtitle = "MSigDb Hallmark")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1), plot.title.position = "plot")
ggsave("MCT.vs.fer.gsea.dotplot.pdf", width = 7, height = 7)
ggsave("MCT.vs.fer.gsea.dotplot.png")

holder <- gsea2plot %>% 
  filter(tests != "MCT.vs.fer" & description != c("HALLMARK_PROTEIN_SECRETION","HALLMARK_GLYCOLYSIS"))

ggplot(data = holder, aes(x=cells, y=pathway))+
  geom_point(aes(size = Gene_ratio, col=FDR))+
  ylab(element_blank())+
  xlab(element_blank())+
  scale_color_viridis_c()+
  ggtitle("Pathways affected in MCT individuals", subtitle = "MSigDb Hallmark")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1), plot.title.position = "plot")
ggsave("control.vs.MCT.gsea.dotplot.pdf", width = 5, height = 7)
ggsave("control.vs.MCT.gsea.dotplot.png")

# Endothelial cells
unwanted <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_MITOTIC_SPINDLE","HALLMARK_DNA_REPAIR","HALLMARK_ADIPOGENESIS")
holder <- gsea2plot %>% 
  filter(cells == "Endothelial") %>%
  mutate(keep = case_when(description %in% unwanted ~ "no", 
                          T ~ "yes")) %>%filter(keep == "yes")
ggplot(data = holder, aes(x=tests, y=pathway))+
  geom_point(aes(size = Gene_ratio, col=FDR))+
  ylab(element_blank())+
  xlab(element_blank())+
  scale_color_viridis_c()+
  ggtitle("Pathways Affected in Endothelial Cells", subtitle = "MSigDb Hallmark")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1), plot.title.position = "plot")
ggsave("Endothelial.gsea.dotplot.pdf", width = 4, height = 5)
ggsave("Endothelial.gsea.dotplot.png")

# All cells
unwanted <- c("HALLMARK_PROTEIN_SECRETION")
holder <- gsea2plot %>% 
  filter(cells == "All_cells") %>%
  mutate(keep = case_when(description %in% unwanted ~ "no", 
                          T ~ "yes")) %>%filter(keep == "yes")
ggplot(data = holder, aes(x=tests, y=pathway))+
  geom_point(aes(size = Gene_ratio, col=FDR))+
  ylab(element_blank())+
  xlab(element_blank())+
  scale_color_viridis_c()+
  ggtitle("Pathways Affected in All Cells", subtitle = "MSigDb Hallmark")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1), plot.title.position = "plot")
ggsave("All_cells.gsea.dotplot.pdf", width = 5, height = 7)
ggsave("All_cells.gsea.dotplot.png")

save.image("toast.Rdata")