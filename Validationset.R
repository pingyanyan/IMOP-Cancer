
tracerxdif<-function(comp){
  install.packages("fst")
  library("fst")
  library("readxl")
  setwd("C://comp_LUAD//validation/TRACERx")
  TPMexp<-read_fst(
    "C://comp_LUAD//validation//TRACERx//Genomic–transcriptomic evolution in lung cancer and metastasis//2022-10-17_rsem_tpm_mat.fst",
    columns = NULL,
    from = 1,
    to = NULL,
    as.data.table = FALSE,
    old_format = FALSE
  )
  mut<-read_fst(
    "C://comp_LUAD//validation//Genomic–transcriptomic evolution in lung cancer and metastasis//20221109_TRACERx421_mutation_table.fst",
    columns = NULL,
    from = 1,
    to = NULL,
    as.data.table = FALSE,
    old_format = FALSE
  )
  cnv<-read_fst(
    "C://comp_LUAD//validation//TRACERx//Evolutionary characterization of lung adenocarcinoma morphology in TRACERx//20221108_sequenza_seg_df_cruk.fst",
    columns = NULL,
    from = 1,
    to = NULL,
    as.data.table = FALSE,
    old_format = FALSE
  )
 clinical<- readRDS("C://comp_LUAD//validation//TRACERx//Evolutionary characterization of lung adenocarcinoma morphology in TRACERx//20221109_TRACERx421_all_patient_df.rds")
  adlung<-clinical[which(clinical$histology_lesion1=="Adenosquamous carcinoma"|clinical$histology_lesion1=="Invasive adenocarcinoma"),1]
  
  #comp<-"CSMD3_PTPRD"
  gene1<-strsplit(comp,"_")[[1]][[1]]
  gene2<-strsplit(comp,"_")[[1]][[2]]
  sampleall<-unique(mut$patient_id)
  sample1<-unique(mut[which(mut$Hugo_Symbol==gene1),2])
  sample2<-unique(mut[which(mut$Hugo_Symbol==gene2),2])
  sample<-intersect(sample1,sample2)
  samplegroup<-NA
  i=1
  for (i in 1:length(sample)){
    ccf1<- mut[which(mut$Hugo_Symbol==gene1&mut$patient_id==sample[i]&mut$func=="exonic"),25]
    ccf2<- mut[which(mut$Hugo_Symbol==gene2&mut$patient_id==sample[i]&mut$func=="exonic"),25]
    if(length(ccf1)>1){ 
      calculate_sum_ccf <- function(mutation_string) {  
        parts <- strsplit(mutation_string, ";")[[1]]   
        ccfs <- as.numeric(sapply(parts, function(part) strsplit(part, ":")[[1]][2]))  
        sum(ccfs)  
      }  
      sum_ccfs <- sapply(ccf1, calculate_sum_ccf)  
      max_sum_index <- which.max(sum_ccfs) 
      ccf1<-ccf1[max_sum_index]
    }
    if(length(ccf2)>1){ 
      calculate_sum_ccf <- function(mutation_string) {  
        parts <- strsplit(mutation_string, ";")[[1]]   
        ccfs <- as.numeric(sapply(parts, function(part) strsplit(part, ":")[[1]][2]))  
        sum(ccfs)  
      }  
      sum_ccfs <- sapply(ccf2, calculate_sum_ccf)  
      max_sum_index <- which.max(sum_ccfs) 
      ccf2<-ccf2[max_sum_index]
    }
    if(length(ccf1)==1&length(ccf2)==1){
      if(is.na(ccf1)==F&is.na(ccf2)==F){
        ccf1<- unlist(strsplit(ccf1,";"))
        ccf2<- unlist(strsplit(ccf2,";"))
        ccf1<-as.character(strsplit2(ccf1, ":")[,2])
        ccf2<-as.character(strsplit2(ccf2, ":")[,2])
        ccf1<-as.numeric(ccf1)
        ccf2<-as.numeric(ccf2)
        region1<-which(ccf1>0)
        region2<-which(ccf2>0)
        if(length(region1)==length(ccf1)&length(region2)==length(ccf2)){
          if(sum(as.numeric(ccf1))>sum(as.numeric(ccf2))){
            sample_name <-"sampleAB"
          }else if(sum(as.numeric(ccf1))<sum(as.numeric(ccf2))){
            sample_name <-"sampleBA"
          }else{
            sample_name <-"sampleAB_BA_unkown"
          }
        }else if(length(region1)==length(ccf1)&length(region2)<length(ccf2)){
          sample_name <-"sampleAB"
        }else if(length(region1)<length(ccf1)&length(region2)==length(ccf2)){
          sample_name <-"sampleBA"
        }else{
          if(length(region1)>length(region2)){
            sample_name <-"sampleAB"
          }else if(length(region1)<length(region2)){
            sample_name <-"sampleBA"
          }else{
            region<-intersect(region1,region2)
            if(length(region)>1){
              if(sum(as.numeric(ccf1[region]))>sum(as.numeric(ccf2[region]))){
                sample_name <-"sampleAB"
              }else if(sum(as.numeric(ccf1[region]))<sum(as.numeric(ccf2[region]))){
                sample_name <-"sampleBA"
              }else{
                sample_name <-"sampleAB_BA_unkown"
              }
            }else if(length(region)==1){
              if(ccf1[region]>ccf2[region]){
                sample_name <-"sampleAB"
              }else if(ccf1[region]<ccf2[region]){
                sample_name <-"sampleBA"
              }else{
                sample_name <-"sampleAB_BA_unkown"
              }
            }else{
              sample_name <-"sampleAB_BA_unkown"
            }
          }
        }
      }else{
        sample_name <-"sampleAB_BA_unkown"
      }
    }else{
      sample_name <-"sampleAB_BA_unkown"}
    sampleid<-cbind(sample[i],sample_name)
    samplegroup<-rbind(samplegroup,sampleid)
  }
  samplegroup<-as.data.frame(samplegroup)
  samplegroup<-samplegroup[-1,]
  sample0A<-setdiff(sample1,sample)
  sample0B<-setdiff(sample2,sample)
  sampleAB<-samplegroup[which(samplegroup[,2]=="sampleAB"),1]
  sampleBA<-samplegroup[which(samplegroup[,2]=="sampleBA"),1]
  sampleunkown<-samplegroup[which(samplegroup[,2]=="sampleAB_BA_unkown"),1]
  sample00<-setdiff(sampleall,c(sample0A,sample0B,sampleAB,sampleBA,sampleunkown)) 

  sample00<-intersect(sample00,adlung)
  sample0A<-intersect(sample0A,adlung)
  sample0B<-intersect(sample0B,adlung)
  sampleAB<-intersect(sampleAB,adlung)
  sampleBA<-intersect(sampleBA,adlung)

  #ssgsea（不合并区域样本）
  TPMexp1<-TPMexp
  rownames(TPMexp1)<-TPMexp1[,1]
  TPMexp1<-TPMexp1[,-1]
  load("C://comp_LUAD//hallmarkall.rdata")
  genesets <- split(msig_t2g_1$gene_symbol, msig_t2g_1$gs_name)
  ssGSEA_hallmark <- gsva(expr = as.matrix(TPMexp1), 
                          gset.idx.list =genesets,
                          method = 'ssgsea',kcdf = 'Gaussian',abs.ranking = TRUE)
  ssgsea_hallmarkall_logTPM<- t(scale(t(ssGSEA_hallmark)))#归一化
  ssgsea_hallmarkall_logTPM[ssgsea_hallmarkall_logTPM< -2] <- -2
  ssgsea_hallmarkall_logTPM[ssgsea_hallmarkall_logTPM>2] <- 2
  normalization <- function(x){return((x-min(x))/(max(x)-min(x)))}  
  ssgsea_hallmarkall_logTPM <- normalization(ssgsea_hallmarkall_logTPM)
  save(ssgsea_hallmarkall_logTPM,file = "C://comp_LUAD//validation//TRACERx//ssgsea_hallmark_logTPM_region.Rdata")
  
  kegg_geneset1 <- qusage::read.gmt("C://comp_LUAD//kegg_hsa.gmt") #返回的是list
  ssGSEA_kegg<- gsva(expr = as.matrix(TPMexp1), 
                     gset.idx.list =kegg_geneset1,
                     method = 'ssgsea',kcdf = 'Gaussian',abs.ranking = TRUE)
  ssgsea_keggall_logTPM<- t(scale(t(ssGSEA_kegg)))#归一化
  ssgsea_keggall_logTPM[ssgsea_keggall_logTPM< -2] <- -2
  ssgsea_keggall_logTPM[ssgsea_keggall_logTPM>2] <- 2
  normalization <- function(x){return((x-min(x))/(max(x)-min(x)))}  
  ssgsea_keggall_logTPM <- normalization(ssgsea_keggall_logTPM)
  save(ssgsea_keggall_logTPM,file = "C://comp_LUAD//validation//TRACERx//ssgsea_kegg_logTPM_region.Rdata")
  
  load("C://comp_LUAD//network_genelist.rdata")
  ssGSEA_network <- gsva(expr = as.matrix(TPMexp1), 
                         gset.idx.list =network,
                         method = 'ssgsea',kcdf = 'Gaussian',abs.ranking = TRUE)
  ssgsea_networkall_logTPM<- t(scale(t(ssGSEA_network)))#归一化
  ssgsea_networkall_logTPM[ssgsea_networkall_logTPM< -2] <- -2
  ssgsea_networkall_logTPM[ssgsea_networkall_logTPM>2] <- 2
  normalization <- function(x){return((x-min(x))/(max(x)-min(x)))}  
  ssgsea_networkall_logTPM <- normalization(ssgsea_networkall_logTPM)
  save(ssgsea_networkall_logTPM,file = "C://comp_LUAD//validation//TRACERx//ssgsea_network_logTPM_region.Rdata")
  
  go_geneset1 <- qusage::read.gmt("C://comp_LUAD//c5.go.v2023.2.Hs.symbols.gmt") 
  ssGSEA_go <- gsva(expr = as.matrix(TPMexp1), 
                    gset.idx.list =go_geneset1,
                    method = 'ssgsea',kcdf = 'Gaussian',abs.ranking = TRUE)
  ssgsea_goall_logTPM<- t(scale(t(ssGSEA_go)))#归一化
  ssgsea_goall_logTPM[ssgsea_goall_logTPM< -2] <- -2
  ssgsea_goall_logTPM[ssgsea_goall_logTPM>2] <- 2
  normalization <- function(x){return((x-min(x))/(max(x)-min(x)))}  
  ssgsea_goall_logTPM <- normalization(ssgsea_goall_logTPM)
  save(ssgsea_goall_logTPM,file ="C://comp_LUAD//validation//TRACERx//ssgsea_go_logTPM_region.Rdata")
  
#多区域样本分组顺序
  a<-colnames(ssgsea_hallmarkall_logTPM)
  sampleAB_1<-a[which(substr(a,1,8)%in%sampleAB)]
  sampleBA_1<-a[which(substr(a,1,8)%in%sampleBA)]
  
  #突变共发生区域
  sampleAB_2<-NA
  i=1
  for(i in 1:length(sampleAB)){
    mut1<-mut[which(mut$patient_id==sampleAB[i]&mut$Hugo_Symbol==gene1),]
    mut2<-mut[which(mut$patient_id==sampleAB[i]&mut$Hugo_Symbol==gene2),]
    ccf1<- mut[which(mut$Hugo_Symbol==gene1&mut$patient_id==sampleAB[i]&mut$func=="exonic"),25]
    ccf2<- mut[which(mut$Hugo_Symbol==gene2&mut$patient_id==sampleAB[i]&mut$func=="exonic"),25]
    if(length(ccf1)>1){ 
      calculate_sum_ccf <- function(mutation_string) {  
        parts <- strsplit(mutation_string, ";")[[1]]   
        ccfs <- as.numeric(sapply(parts, function(part) strsplit(part, ":")[[1]][2]))  
        sum(ccfs)  
      }  
      sum_ccfs <- sapply(ccf1, calculate_sum_ccf)  
      max_sum_index <- which.max(sum_ccfs) 
      ccf1<-ccf1[max_sum_index]
    }
    if(length(ccf2)>1){ 
      calculate_sum_ccf <- function(mutation_string) {  
        parts <- strsplit(mutation_string, ";")[[1]]   
        ccfs <- as.numeric(sapply(parts, function(part) strsplit(part, ":")[[1]][2]))  
        sum(ccfs)  
      }  
      sum_ccfs <- sapply(ccf2, calculate_sum_ccf)  
      max_sum_index <- which.max(sum_ccfs) 
      ccf2<-ccf2[max_sum_index]
    }
    ccf1<- unlist(strsplit(ccf1,";"))
    ccf2<- unlist(strsplit(ccf2,";"))
    ccf_1<-as.character(strsplit2(ccf1, ":")[,2])
    ccf_2<-as.character(strsplit2(ccf2, ":")[,2])
    index<-intersect(which(ccf_1!=0),which(ccf_2!=0))
    region<-as.character(strsplit2(ccf1, ":")[,1])
    region<-gsub("[.]","-",region[index])
    samplename<-paste0(sampleAB[i],"_",region)
    sampleAB_2<-c(sampleAB_2,samplename)
  }
  sampleAB_2<-sampleAB_2[-1]
  
  sampleBA_2<-NA
  i=1
  for(i in 1:length(sampleBA)){
    mut1<-mut[which(mut$patient_id==sampleBA[i]&mut$Hugo_Symbol==gene1),]
    mut2<-mut[which(mut$patient_id==sampleBA[i]&mut$Hugo_Symbol==gene2),]
    ccf1<- mut[which(mut$Hugo_Symbol==gene1&mut$patient_id==sampleBA[i]&mut$func=="exonic"),25]
    ccf2<- mut[which(mut$Hugo_Symbol==gene2&mut$patient_id==sampleBA[i]&mut$func=="exonic"),25]
    if(length(ccf1)>1){ 
      calculate_sum_ccf <- function(mutation_string) {  
        parts <- strsplit(mutation_string, ";")[[1]]   
        ccfs <- as.numeric(sapply(parts, function(part) strsplit(part, ":")[[1]][2]))  
        sum(ccfs)  
      }  
      sum_ccfs <- sapply(ccf1, calculate_sum_ccf)  
      max_sum_index <- which.max(sum_ccfs) 
      ccf1<-ccf1[max_sum_index]
    }
    if(length(ccf2)>1){ 
      calculate_sum_ccf <- function(mutation_string) {  
        parts <- strsplit(mutation_string, ";")[[1]]   
        ccfs <- as.numeric(sapply(parts, function(part) strsplit(part, ":")[[1]][2]))  
        sum(ccfs)  
      }  
      sum_ccfs <- sapply(ccf2, calculate_sum_ccf)  
      max_sum_index <- which.max(sum_ccfs) 
      ccf2<-ccf2[max_sum_index]
    }
    ccf1<- unlist(strsplit(ccf1,";"))
    ccf2<- unlist(strsplit(ccf2,";"))
    ccf_1<-as.character(strsplit2(ccf1, ":")[,2])
    ccf_2<-as.character(strsplit2(ccf2, ":")[,2])
    index<-intersect(which(ccf_1!=0),which(ccf_2!=0))
    region<-as.character(strsplit2(ccf1, ":")[,1])
    region<-gsub("[.]","-",region[index])
    samplename<-paste0(sampleBA[i],"_",region)
    sampleBA_2<-c(sampleBA_2,samplename)
  }
  sampleBA_2<-sampleBA_2[-1]
  save(sampleAB_2,sampleBA_2,file=paste0("C://comp_LUAD//validation//",comp,"//",comp,"_tracerx_sample2.rdata"))
  
}

cell_dif<-function(comp){
  #comp<-"CSMD3_PTPRD"
  gene1<-strsplit(comp,"_")[[1]][[1]]
  gene2<-strsplit(comp,"_")[[1]][[2]]
  setwd("C://comp_LUAD//validation/CELL")
  
  mutall<-read.csv("C://comp_LUAD//validation//CELL//OmicsSomaticMutations.csv")
  cell<-unique(mutall$ModelID)
  id<-read.csv("C://comp_LUAD//validation//Model.csv")
  cell_lung<-id[which(id$OncotreeSubtype=="Lung Adenocarcinoma"),1]
  cell_lung<-intersect(cell_lung,cell)
  mut_lung<-mutall[ which(mutall$ModelID %in% cell_lung),]
  ABSOLUTE<-read.csv("C://comp_LUAD//validation//CELL//ABSOLUTE.csv")
  samplegroup<-NA
  i=1
  for(i in 1:length(cell_lung)){
    c<-subset(mut_lung, ModelID == as.character(cell_lung[i]))
    d<-subset(ABSOLUTE, depMapID == as.character(cell_lung[i]))
    maf<-c[,c(15,1,2,3,4,7,8)]
    cnv<-d[,c(2,3,4,7,8)]
    maf$Chrom =substring(maf$Chrom,4)
    maf$Chrom = as.character(maf$Chrom)
    maf$Pos = as.numeric(maf$Pos)
    cnv$Chromosome = as.character(cnv$Chromosome)
    cnv$Start = as.numeric(cnv$Start)
    cnv$End = as.numeric(cnv$End)
    npos <- length(maf$Pos)
    local_cn_a1<- vector(mode="numeric",length=npos)
    local_cn_a2<- vector(mode="numeric",length=npos)
    Chrs <- unique(maf$Chrom)
    for(chr in Chrs){
      t1.chr.matches <- which(maf$Chrom==chr)
      t2.chr.compares <- which(cnv$Chromosome==chr)
      for(match in t1.chr.matches){
        cn_a1=0
        cn_a2=0
        candidate = as.numeric(as.vector(maf$Pos)[match])
        for(compare in t2.chr.compares){
          comp.start <- as.numeric(as.vector(cnv$Start)[compare])
          comp.stop <- as.numeric(as.vector(cnv$End)[compare])
          if(candidate>=comp.start&&candidate<=comp.stop){
            cn_a1 = as.numeric(as.vector(cnv$Modal_HSCN_1)[compare])
            cn_a2 = as.numeric(as.vector(cnv$Modal_HSCN_2)[compare])
            break
          }
        }
        local_cn_a1[match] = cn_a1
        local_cn_a2[match] = cn_a2
      }
    }
    Tab1 <- as.data.frame(cbind(maf,local_cn_a1,local_cn_a2))
    colnames(Tab1)<-c("Hugo_Symbol","Chromosome","Start_position","Reference_Allele","Tumor_Seq_Allele2","t_ref_count","t_alt_count","local_cn_a1","local_cn_a2")
    outfile<-paste("C://comp_LUAD//validation//CELL//CELL_NDTinput//",cell_lung[i],"_input",".txt",sep="")
    write.table(Tab1,file=outfile,quote = F,sep="\t",row.names = F)
  }
  
  purity<-as.data.frame(cell_lung)
  purity$purity<-1
  write.table(purity,file = "./CELL_NDTinput/purity.txt",quote = FALSE, sep = "\t",row.names = FALSE,col.names=FALSE)
  
  #ascat
  b<-ABSOLUTE[,c(2,3,4,8,7,20)]
  colnames(b)<-c("Chromosome","Start_position","End_Position","A2_CN","A1_CN","sample")
  b$Chromosome = as.character(b$Chromosome)
  b$Start_position = as.numeric(b$Start_position)
  b$End_Position = as.numeric(b$End_Position)
  b$A1_CN = as.numeric(b$A1_CN)
  b$A2_CN = as.numeric(b$A2_CN)
  b$sample = as.character(b$sample)
  i=1
  for(i in 1:nrow(purity)){
    d<-subset(b, sample == as.character(purity[i,1]))
    d<-d[,c(1,2,3,5,4)]
    outfile<-paste("C://comp_LUAD//validation//CELL//CELL_NDTinput//",purity[i,1],"_ascat",".txt",sep="")
    write.table(d,file=outfile,quote = F,sep="\t",row.names = F)
  }
  
  #每个肺腺癌细胞系运行PhylogicNDT
  
  #samplegroup
  samplegroup<-NA
  gene1<-"CSMD3"
  gene2<-"PTPRD"
  i=1
  for(i in 1:length(cell_lung)){
    mut<-mut_lung[which(mut_lung$ModelID==cell_lung[i]),]
    geneall<-unique(mut$HugoSymbol)
    if (gene1 %in% geneall & gene2 %in% geneall) {
      sample_name <- "sampleAB_BA"
    } else if (gene1 %in% geneall) {
      sample_name <- "sample0A"
    } else if (gene2 %in% geneall) {
      sample_name <- "sample0B"
    } else {
      sample_name <- "sample_no"
    }
    samplegroup<-c(samplegroup,sample_name)
  }
  samplegroup<-samplegroup[-1]
  samplegroup<-rbind(cell_lung,samplegroup)
  samplegroup<-as.data.frame(t(samplegroup))
  cell_00<-samplegroup[which(samplegroup$samplegroup=="sample_no"),1]
  cell_0A<-samplegroup[which(samplegroup$samplegroup=="sample0A"),1]
  cell_0B<-samplegroup[which(samplegroup$samplegroup=="sample0B"),1]
  cell_AB_BA<-samplegroup[which(samplegroup$samplegroup=="sampleAB_BA"),1]
  
  i=1
  cell_AB<-NA
  cell_BA<-NA
  cell_unkown<-NA
  for(i in 1:length(cell_AB_BA)){
    a<-read.table(paste0("C://comp_LUAD//validation//CELL//NDT//",cell_AB_BA[i],"//",cell_AB_BA[i],".mut_ccfs.txt"),
                  header=T,sep="\t",quote="")
    if(max(as.numeric(a[which(a$Hugo_Symbol==gene1),17]))>max(as.numeric(a[which(a$Hugo_Symbol==gene2),17]))){
      cell_AB<-c(cell_AB,cell_AB_BA[i])
    }else if(max(as.numeric(a[which(a$Hugo_Symbol==gene1),17]))<max(as.numeric(a[which(a$Hugo_Symbol==gene2),17]))){
      cell_BA<-c(cell_BA,cell_AB_BA[i])
    }else{
      cell_unkown<-c(cell_unkown,cell_AB_BA[i])
    }
  }
  cell_AB<-cell_AB[-1]
  cell_BA<-cell_BA[-1]
  save(cell_00,cell_0A,cell_0B,cell_AB,cell_BA,file = "cell_A_B.rdata")
}

oncosg_dif<-function(comp){
  #comp<-"CSMD3_PTPRD"
  gene1<-strsplit(comp,"_")[[1]][[1]]
  gene2<-strsplit(comp,"_")[[1]][[2]]
  setwd("C://comp_LUAD//validation/ONCOSG")
  maf<-read.table("C://comp_LUAD//validation//ONCOSG//oncosg_mutations.txt",header=T,sep="\t",quote="")
  seg<-read.table("C://comp_LUAD//validation//ONCOSG//oncosg//oncosg_cna_hg19.seg",header=T,sep = "\t",quote="")
  colnames(seg)=c("Sample","Chromosome","Start","End","Num_Probes","Segment_Mean")
  library(numDeriv)
  library(DoAbsolute)
  library(ABSOLUTE)
  result<-DoAbsolute(Seg = seg,
                     Maf =  maf,
                     platform = "SNP_6.0",
                     copy.num.type = "total",
                     results.dir = "test_oncosg",
                     nThread = 2,
                     keepAllResult = TRUE,
                     verbose = TRUE)
  
  i=1
  oncosg_AB<-NA
  oncosg_BA<-NA
  oncosg_same<-NA
  filenames <- list.files("C://comp_LUAD//validation//ONCOSG//test_oncosg//maf", full.names = TRUE)
  sample<-substr(filenames,46,nchar(filenames)-12)
  i=1
  for(i in 1:length(filenames)){
    ccf<-read.table(filenames[i],header=T,sep = "\t",quote="")
    a<-ccf[which(ccf$Hugo_Symbol==gene1),31]
    b<-ccf[which(ccf$Hugo_Symbol==gene2),31]
    if(length(a)>0&length(b)>0){
      if(max(a)>max(b)){
        oncosg_AB<-c(oncosg_AB,sample[i])
      }else if(max(a)<max(b)){
        oncosg_BA<-c(oncosg_BA,sample[i])
      }else{
        a1<-ccf[which(ccf$Hugo_Symbol==gene1),30] 
        b1<-ccf[which(ccf$Hugo_Symbol==gene2),30]
        if(max(a1)>max(b1)){
          oncosg_AB<-c(oncosg_AB,sample[i])
        }else if(max(a1)<max(b1)){
          oncosg_BA<-c(oncosg_BA,sample[i])
        }else{
          oncosg_same<-c(oncosg_same,sample[i])
        }
      }
    }
  }
  oncosg_AB<-oncosg_AB[-1]
  oncosg_BA<-oncosg_BA[-1]
  save(oncosg_AB,oncosg_BA,file="oncosg_AB_BA.rdata")

}

su2clcdif<-function(comp){
  #comp<-"CSMD3_PTPRD"
  gene1<-strsplit(comp,"_")[[1]][[1]]
  gene2<-strsplit(comp,"_")[[1]][[2]]
  library(readxl)
  setwd("C://comp_LUAD//validation/SU2CLC")
  su2clc_clinical<-read_excel("C://comp_LUAD//validation/SU2CLC//41588_2023_1355_MOESM3_ESM.xlsx", sheet = "Table_S1_Clinical_Annotations",skip=2)
  su2clc_mut<-read_excel("C://comp_LUAD//validation/SU2CLC//41588_2023_1355_MOESM3_ESM.xlsx", sheet = "Table_S2_Absolute_Maf")
  #顺序样本
  sample<-as.character(t(su2clc_clinical[which(su2clc_clinical$Histology_Harmonized=="Adeno"),8]))
  sample1<-intersect(as.character(t(su2clc_mut[which(su2clc_mut$Hugo_Symbol==gene1),16])),
                     as.character(t(su2clc_mut[which(su2clc_mut$Hugo_Symbol==gene2),16])))
  su2clc_00<-NA
  su2clc_0A<-NA
  su2clc_0B<-NA
  su2clc_AB<-NA
  su2clc_BA<-NA
  su2clc_unknown<-NA
  sample<-paste0(sample,"-T1")
  i=1
  for (i in 1:length(sample)){
    mut<-su2clc_mut[which(su2clc_mut$Tumor_Sample_Barcode ==sample[i]),]
    if(nrow(mut)>0){
      geneall<-unique(mut$Hugo_Symbol)
      if(gene1%in%geneall&gene2%in%geneall){
        ccf1<-max(mut[which(mut$Hugo_Symbol==gene1),52])
        ccf2<-max(mut[which(mut$Hugo_Symbol==gene2),52])
        if(length(na.omit(ccf1))>0&length(na.omit(ccf2))>0){
          if(ccf1>ccf2){
            su2clc_AB<-c(su2clc_AB,sample[i])
          }else if(ccf1<ccf2){
            su2clc_BA<-c(su2clc_BA,sample[i])
          }else{
            ccf111<-max(mut[which(mut$Hugo_Symbol==gene1),34])
            ccf222<-max(mut[which(mut$Hugo_Symbol==gene2),34])
            if(ccf111>ccf222){
              su2clc_AB<-c(su2clc_AB,sample[i])
            }else if(ccf111<ccf222){
              su2clc_BA<-c(su2clc_BA,sample[i])
            }else{
              su2clc_unknown<-c(su2clc_unknown,sample[i])
            }
          }
        }
      }
      else if(gene1%in%geneall){
        su2clc_0A<-c(su2clc_0A,sample[i])
      }else if(gene2%in%geneall){
        su2clc_0B<-c(su2clc_0A,sample[i])
      }else{
        su2clc_00<-c(su2clc_00,sample[i])
      }
    }
  }
  su2clc_00<-su2clc_00[-1]
  su2clc_0A<-su2clc_0A[-1]
  su2clc_0B<-su2clc_0B[-1]
  su2clc_AB<-su2clc_AB[-1]
  su2clc_BA<-su2clc_BA[-1]
  save(su2clc_00,su2clc_0A,su2clc_0B,su2clc_AB,su2clc_BA,file="su2clc_AB_BA.rdata")
}
