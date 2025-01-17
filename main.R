#main
install.packages("BiocManager")
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
BiocManager::install('UCSCXenaTools',force = TRUE)

library("TCGAbiolinks")
library("dplyr")
library("limma")
library("biomaRt")
library("SummarizedExperiment")
library("stringr")
library("ggplotify")
library("patchwork")
library("cowplot")
library("dplyr")
library('edgeR')
library('tidyr')
library("survminer")
library("ggplot2")
library("survival")
library("msigdbr")
library("clusterProfiler")
library(reshape2)
library(ggpubr)
library(org.Hs.eg.db)
library(GSVA)
library(RColorBrewer)
library(ggpmisc)
library(ggsci)
library(data.table)
library(tidyverse)
library(Biobase)
library(pheatmap)
library(UCSCXenaTools)

##################cancer types in TCGA and their corresponding normal tissues in GTEX
cancer<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
          "KICH","KIRC","KIRP","LAML","LGG","LIHC","LUSC","MESO","OV","PAAD",
          "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC",
          "UCS","UVM")
cancer_tissue<-c( "Adrenal Gland","Bladder","Breast","Cervix Uteri","Liver","Colon","Blood","Esophagus","Brain","Skin",
                 "Kidney","Kidney","Kidney","Blood" ,"Brain","Liver","Lung", "Breast", "Ovary","Pancreas",
                 "Adrenal Gland","Prostate","Colon","Muscle","Skin","Stomach","Testis","Thyroid","Breast" ,"Uterus" ,
                 "Uterus" ,"Liver")#"UVM" has no corresponding ocular tissue and is replaced by liver tissue from the metastatic site

##################Download TCGA mutation spectrum copy number spectrum for each cancer and integrate it into PhylogicNDT input file
source("NDTinput.R")
i=1
for (i in 1:length(cancer)){
  #NDTinput
  get_maf(cancer[i])
  get_clinical(cancer[i])
  sampleinput(cancer[i])
}

##################Download and install PhylogicNDT, run it in the Anaconda environment
#conda activate python27
#d:
#cd PhylogicNDT
#python p.py
#python ptree.py
#python ptime.py


##################Screening for key mutation order pairs
source("get_driver.R")
i=1
for (i in 1:length(cancer)){
  #samplegroup
  #Select mutations with a mutation frequency greater than 0.05 and all copy number variations, and combine them in pairs
  #Calculate the number of grouped samples for each mutation pair separately
  ##sample where neither mutation exists：sample00
  ##sample with only one mutation：sampleA,sampleB
  ##Two mutations occur, but in a different order of occurrence：sampleAB,sampleBA
  samplegroup_mut_mut(cancer[i])
  samplegroup_mut_scan(cancer[i])
  samplegroup_scan_scan(cancer[i])
  #ssGSEA
  get_gtex(cancer_tissue[i])
  get_cancer_fpkm(cancer[i])
  ssgsea_hallmark(cancer[i])
  ssgsea_kegg(cancer[i])
  ssgsea_K_network(cancer[i])
  ssgsea_GO(cancer[i])
  #Screening for key mutation order pairs
  get_driver(cancer[i])
}

##################Identify different occurrence order samples of mutation pairs in the validation set(LUAD)
comp<-""
source("Validationset.rdata")
tracerxdif(comp)
cell_dif(comp)
oncosg_dif(comp)
su2clcdif(comp)