####Header############################################################
# Project: Prins 2023: Identification of differential gene expression 
# responses in pulmonary cells from bulk RNA-Seq using deconvolution
# Code goal: deconvolution of bulk RNAseq libraries
# RI Analyst: Natalia Calixto Mancipe

## Libraries----------------------------------------------------------
.libPaths(c("/panfs/jay/groups/0/riss/calix015/myR_lib/",.libPaths()))
library(tidyverse)
library(DWLS)
library(kableExtra)

## Setup--------------------------------------------------------------
# This session was run in OpenOnDemand using 2 cores and 1 node, R/4.1.0
setwd("~/2023/prins/")
indir="~/2023/prins/bulk_data/"
outdir = "~/2023/prins/deconv"
if(!file.exists(outdir)) dir.create(outdir, recursive = TRUE)

#load("scPrep.RData")
load("deconv.Rdata")

## Bulk RNA-Seq --------------------
# Load bulkRNAseq data
all_reads <- read_delim(paste0(indir, "subread_counts_gene_symbol.txt"), delim = "\t", col_names = T, comment = "#")

# Separate control, MCT libraries and gene names
# Put gene names in the first column 
colnames(all_reads)
control.libs <- all_reads[,c(1,which(str_detect(colnames(all_reads), pattern = "control")))]
mct.libs <- all_reads[,c(1,6:length(colnames(all_reads)))]

## Deconvolve controls ------
load("sc_data/Control_MAST/Sig.RData")
allCounts_DWLS<-NULL
counts <- data.matrix(control.libs[2:dim(control.libs)[2]])
rownames(counts) <- control.libs$GeneName

for(j in 1:(dim(counts)[2])){
  S<-Sig
  Bulk<-counts[,j]
  Genes<-intersect(rownames(S),rownames(counts))
  B<-Bulk[Genes]
  S<-S[Genes,]
  sample <- colnames(counts)[j]
  solDWLS<-solveDampenedWLS(S,B)
  file.name <- paste0(outdir,"/DWLS_",sample,".txt")
  write.table(solDWLS, file = file.name, append = F, quote = F, sep = "\t", col.names = F)
}

## Deconvolve MCTs with MCT sc. reference ------
load("sc_data/MCT_MAST/Sig.RData")
allCounts_DWLS<-NULL
counts <- data.matrix(mct.libs[2:dim(mct.libs)[2]])
rownames(counts) <- mct.libs$GeneName

for(j in 1:(dim(counts)[2])){
  S<-Sig
  Bulk<-counts[,j]
  Genes<-intersect(rownames(S),rownames(counts))
  B<-Bulk[Genes]
  S<-S[Genes,]
  sample <- colnames(counts)[j]
  solDWLS<-solveDampenedWLS(S,B)
  file.name <- paste0(outdir,"/DWLS_",sample,".txt")
  write.table(solDWLS, file = file.name, append = F, quote = F, sep = "\t", col.names = F)
}

## Deconvolve MCTs with control sc. reference ------
load("sc_data/Control_MAST/Sig.RData")
allCounts_DWLS<-NULL
counts <- data.matrix(mct.libs[2:dim(mct.libs)[2]])
rownames(counts) <- mct.libs$GeneName

for(j in 1:(dim(counts)[2])){
  S<-Sig
  Bulk<-counts[,j]
  Genes<-intersect(rownames(S),rownames(counts))
  B<-Bulk[Genes]
  S<-S[Genes,]
  sample <- colnames(counts)[j]
  solDWLS<-solveDampenedWLS(S,B)
  file.name <- paste0(outdir,"/DWLS_scControl_",sample,".txt")
  write.table(solDWLS, file = file.name, append = F, quote = F, sep = "\t", col.names = F)
}

## Figures ----
# Proportion data and make df, then numeric matrix
file_names <- list.files(outdir, pattern = "^DWLS_NV_", full.names = T)

Prop <- NULL

for (i in 1:length(file_names)){
  x.file <- file_names[i]
  test <- str_split(x.file, pattern = "/", simplify = T)
  test <- test[,dim(test)[2]]
  test <- str_replace(test, pattern = ".txt",replacement = "")
  test <- str_replace(test, pattern = "DWLS_NV_",replacement = "")
  
  a <- read_delim(x.file, delim = "\t", col_names = F) 
  colnames(a) <- c("cell_type",test)
  if(i == 1){
    Prop <- a
  }
  if(i != 1){
    Prop <- cbind(Prop, a[,2])
  }
}

# "a" will have all the deconvolution results from all files
a <- Prop %>% pivot_longer(-cell_type, names_to = "sample",values_to = "prop")

a <- a %>%
  mutate(treatment = case_when(str_detect(sample, pattern = "control") ~ "control",
                               str_detect(sample, pattern = "fer") ~ "ferrostatin",
                               T ~ "MCT"),
         pc.cell = prop*100)

# Barplot
ggplot(a, aes(x=sample,y=pc.cell))+
  geom_col(aes(group=cell_type, fill=cell_type), col="grey90")+
  ylab("Cell type (%)")+
  xlab(element_blank())+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(outdir,"/deconv_barplot.pdf"))

# Boxplot for each cell type
c.type <- unique(c.types)
for(i in c.type){
  b <- a%>%filter(cell_type == i)
  plot2save <- ggplot(b, aes(x=treatment, y=pc.cell))+
    geom_boxplot()+
    geom_point()+
    ylab("Cell type abundance (%)")+
    ggtitle(i)+
    theme_bw()+
    theme(axis.title.x = element_blank())
    
  ggsave(filename = paste0(outdir,"/",i,"deconv_boxplot.pdf"),
        plot = plot2save, width = 4, height = 4)
}

# Deconvolution results table
Prop <- a %>%
  group_by(pick(treatment,cell_type)) %>%
  summarise(avg = mean(pc.cell),
            sd = sd(pc.cell))

Prop <- Prop %>%
  pivot_wider(names_from = treatment, values_from = c(avg, sd))

fancy.table <- Prop %>%
  arrange(desc(avg_control))%>%
  select(cell_type,avg_control, sd_control, avg_MCT,sd_MCT, avg_ferrostatin, sd_ferrostatin) %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = F)
save_kable(fancy.table, "deconv/table.pdf")

save.image("deconv.Rdata")
