#samplegroup
samplegroup_mut_mut<-function(cancer){
  #mutations with a mutation frequency greater than 0.05
  library(maftools)
  load(paste0("E://comp_allcancer//comp_",cancer,"//TCGA-",cancer,"_SNP.rdata"))
  data<-maf@data
  freq<-NA
  i=1
  mutall<-unique(data$Hugo_Symbol)
  sampleall<-unique(data$Tumor_Sample_Barcode)
  for(i in 1:length(mutall)){
    sample<-data[which(data$Hugo_Symbol==mutall[i]),16]
    mutfreq<-length(unique(sample$Tumor_Sample_Barcode))/length(sampleall)
    freq1<-cbind(mutall[i],mutfreq)
    freq<-rbind(freq,freq1)
  }
  freq<-freq[-1,]
  freq<-as.data.frame(freq)
  mut<-as.character(freq[which(freq$mutfreq>0.05),1])
  mut<-mut[which(mut!="TTN")]
  save(mut,file=paste0("C://comp_allcancer//",cancer,"//mut_freq_0.05.rdata"))
  
  file1=paste0("E://comp_allcancer//comp_",cancer,"//NDTinput//purity_",cancer,".txt")
  purity<-read.table(file1,header=T,sep = "\t")
  ccfall<-NA
  j=1
  for(j in 1:nrow(purity)){
    fileall<-list.files(paste0("E://comp_allcancer//comp_",cancer,"//NDT//",purity[j,1]), full.names = TRUE)
    file3<-paste0("E://comp_allcancer//comp_",cancer,"//NDT//",purity[j,1],"/",purity[j,1],".mut_ccfs.txt")
    if(length(which(fileall%in%file3))>0){
      ccf<-read.table(file3,header=T,sep = "\t")
      ccfall<-rbind(ccfall,ccf)
    }
  }
  ccfall<-ccfall[-1,]
  ccfall$Hugo_Symbol<-as.character(strsplit2(ccfall$Hugo_Symbol," _ ")[,1])
  save(ccfall,file="C://comp_allcancer//",cancer,"//ccfall.rdata")
  
  treeall<-NA
  j=1
  for(j in 1:nrow(purity)){
    fileall<-list.files(paste0("E://comp_allcancer//comp_",cancer,"//NDT//",purity[j,1]), full.names = TRUE)
    file5<-paste0("E://comp_allcancer//comp_",cancer,"//NDT//",purity[j,1],"/",purity[j,1],"_build_tree_posteriors.tsv")
    if(length(which(fileall%in%file5))>0){
      tree<-read.table(file5,header=T,sep = "\t")
      tree$sample<-purity[j,1]
      if(nrow(tree)>0){
        treeall<-rbind(treeall,tree)
      }
    }
  }
  treeall<-treeall[-1,]
  save(treeall,file="C://comp_allcancer//",cancer,"//treeall.rdata")
  
  timeall<-NA
  j=1
  for(j in 1:nrow(purity)){
    fileall<-list.files(paste0("E://comp_allcancer//comp_",cancer,"//NDT//",purity[j,1]), full.names = TRUE)
    file5<-paste0("E://comp_allcancer//comp_",cancer,"//NDT//",purity[j,1],"/",purity[j,1],".timing.tsv")
    if(length(which(fileall%in%file5))>0){
      tim<-read.table(file5,header=T,sep = "\t")
      if(nrow(tim)>0){
        timeall<-rbind(timeall,tim)
      }
    }
  }
  timeall<-timeall[-1,]
  save(treeall,file="C://comp_allcancer//",cancer,"//timeall.rdata")
  
  setwd(paste0("C://comp_allcancer//",cancer,"//samplegroup"))
  
  i=1
  for(i in 1:length(mut)){
    gene1<-mut[i]
    m=sum(i,1)
    for (m in sum(i,1):length(mut)){
      gene2<-mut[m]
      sample1<-ccfall[which(ccfall$Hugo_Symbol==gene1),1]
      sample2<-ccfall[which(ccfall$Hugo_Symbol==gene2),1]
      sample<-intersect(sample1,sample2)
      if(length(sample)>30){
        sampleAB<-NA
        sampleBA<-NA
        sampleunkown<-NA
        j=1
        for(j in 1:length(sample)){
          ccf<-ccfall[which(ccfall$Patient_ID==sample[j]),]
          ccf_1<-max(ccf[which(ccf$Hugo_Symbol==gene1),17])
          ccf_2<-max(ccf[which(ccf$Hugo_Symbol==gene2),17])
          tree<-treeall[which(treeall$sample==sample[j]),]
          tree<-strsplit(tree[1,3],",")[[1]]
          if(length(tree)>1){
            tree<-tree[-length(tree)]
          }else{
            tree<-tree
          }
          cluster_1<-ccf[which(ccf$Hugo_Symbol==gene1),14][which.max(ccf[which(ccf$Hugo_Symbol==gene1),17])]
          cluster_2<-ccf[which(ccf$Hugo_Symbol==gene2),14][which.max(ccf[which(ccf$Hugo_Symbol==gene2),17])]
          if(length(tree)>0){
            branch_points <- strsplit(tree, "-")
            start_points <- sapply(branch_points, function(x) as.numeric(x[1]))
            end_points <- sapply(branch_points, function(x) as.numeric(x[2]))
            if (cluster_1%in% start_points & cluster_2 %in% end_points) {
              sample_name <-"sampleAB"
              sampleAB<-c(sampleAB,sample[j])
            }else if (cluster_2%in% start_points & cluster_1 %in% end_points){
              sample_name <-"sampleBA"
              sampleBA<-c(sampleBA,sample[j])
            }else if(cluster_1==cluster_2&ccf_1>ccf_2){
              sample_name <-"sampleAB"
              sampleAB<-c(sampleAB,sample[j])
            }else if(cluster_1==cluster_2&ccf_1<ccf_2){
              sample_name <-"sampleBA"
              sampleBA<-c(sampleBA,sample[j])
            }else{
              sample_name <-"sampleAB_BA_unkown"
              sampleunkown<-c(sampleunkown,sample[j])
            }
          }
        }
        sampleAB<-sampleAB[-1]
        sampleBA<-sampleBA[-1]
        sampleunkown<-sampleunkown[-1]
        
        sample1<-timeall[which(timeall$Event.Name==gene1),1]
        sample2<-timeall[which(timeall$Event.Name==gene2),1]
        sampleall<-unique(timeall$Patient_ID)
        sample0A<-setdiff(sample1,c(sampleAB,sampleBA))
        sample0B<-setdiff(sample2,c(sampleAB,sampleBA))
        sample00<-setdiff(sampleall,c(sample0A,sample0B,sampleAB,sampleBA))
        
        if(length(sampleAB)>15&length(sampleBA)>15){
          file7=paste0(gene1,"_",gene2,"_sample_A_B.rdata")
          save(sample00,sample0A,sample0B,sampleAB,sampleBA,sampleunkown,file=file7)
        }
      }
    }
  }
  
}

samplegroup_mut_scan<-function(cancer){
  ssetwd(paste0("C://comp_allcancer//",cancer,"//samplegroup"))
  
  load("mut_freq_0.05.rdata")
  
  scan<-c("WGD","gain_1p","loss_1p","gain_1q","loss_1q",
             "gain_2p","loss_2p","gain_2q","loss_2q",
             "gain_3p","loss_3p","gain_3q","loss_3q",
             "gain_4p","loss_4p","gain_4q","loss_4q",
             "gain_5p","loss_5p","gain_5q","loss_5q",
             "gain_6p","loss_6p","gain_6q","loss_6q",
             "gain_7p","loss_7p","gain_7q","loss_7q",
             "gain_8p","loss_8p","gain_8q","loss_8q",
             "gain_9p","loss_9p","gain_9q","loss_9q",
             "gain_10p","loss_10p","gain_10q","loss_10q",
             "gain_11p","loss_11p","gain_11q","loss_11q",
             "gain_12p","loss_12p","gain_12q","loss_12q",
             "gain_13p","loss_13p","gain_13q","loss_13q",
             "gain_14p","loss_14p","gain_14q","loss_14q",
             "gain_15p","loss_15p","gain_15q","loss_15q",
             "gain_16p","loss_16p","gain_16q","loss_16q",
             "gain_17p","loss_17p","gain_17q","loss_17q",
             "gain_18p","loss_18p","gain_18q","loss_18q",
             "gain_19p","loss_19p","gain_19q","loss_19q",
             "gain_20p","loss_20p","gain_20q","loss_20q",
             "gain_21p","loss_21p","gain_21q","loss_21q",
             "gain_22p","loss_22p","gain_22q","loss_22q")
  
  file1=paste0("E://comp_allcancer//comp_",cancer,"//NDTinput//purity_",cancer,".txt")
  purity<-read.table(file1,header=T,sep = "\t")
  load("C://comp_allcancer//",cancer,"//timeall.rdata")
  
  setwd(paste0("C://comp_allcancer//",cancer,"//samplegroup"))
  i=1
  for(i in 1:length(mut)){
    gene1<-mut[i]
    m=1
    for (m in 1:length(scan)){
      gene2<-scan[m]
      sample1<-timeall[which(timeall$Event.Name==gene1),1]
      sample2<-timeall[which(timeall$Event.Name==gene2),1]
      sample<-intersect(sample1,sample2)
      if(length(sample)>30){
        sampleAB<-NA
        sampleBA<-NA
        sampleunkown<-NA
        j=1
        for(j in 1:length(sample)){
          tim<-timeall[timeall$Patient_ID==sample[j],]
          tim_1<-min(tim[which(tim$Event.Name==gene1),13])
          tim_2<-min(tim[which(tim$Event.Name==gene2),13])
          #Sometimes PhylogicNDT cannot infer an accurate false time probability density curve,result a mean of 0.05,removing these interference values
          if (tim_1!=0.05&tim_2!=0.05&tim_1<tim_2) {
            sample_name <-"sampleAB"
            sampleAB<-c(sampleAB,sample[j])
          }else if (tim_1!=0.05&tim_2!=0.05&tim_1>tim_2){
            sample_name <-"sampleBA"
            sampleBA<-c(sampleBA,sample[j])
          }else{
            sample_name <-"sampleAB_BA_unkown"
            sampleunkown<-c(sampleunkown,sample[j])
          }
        }
        sampleAB<-sampleAB[-1]
        sampleBA<-sampleBA[-1]
        sampleunkown<-sampleunkown[-1]
        
        sample1<-timeall[which(timeall$Event.Name==gene1),1]
        sample2<-timeall[which(timeall$Event.Name==gene2),1]
        sampleall<-unique(timeall$Patient_ID)
        sample0A<-setdiff(sample1,c(sampleAB,sampleBA))
        sample0B<-setdiff(sample2,c(sampleAB,sampleBA))
        sample00<-setdiff(sampleall,c(sample0A,sample0B,sampleAB,sampleBA))
        
        if(length(sampleAB)>15&length(sampleBA)>15){
          file7=paste0(gene1,"_",gene2,"_sample_A_B.rdata")
          save(sampleAB,sampleBA,sampleunkown,file=file7)
        }
      }
    }
  }
}      

samplegroup_scan_scan<-function(cancer){
  setwd(paste0("C://comp_allcancer//",cancer,"//samplegroup"))
  scan<-c("WGD","gain_1p","loss_1p","gain_1q","loss_1q",
          "gain_2p","loss_2p","gain_2q","loss_2q",
          "gain_3p","loss_3p","gain_3q","loss_3q",
          "gain_4p","loss_4p","gain_4q","loss_4q",
          "gain_5p","loss_5p","gain_5q","loss_5q",
          "gain_6p","loss_6p","gain_6q","loss_6q",
          "gain_7p","loss_7p","gain_7q","loss_7q",
          "gain_8p","loss_8p","gain_8q","loss_8q",
          "gain_9p","loss_9p","gain_9q","loss_9q",
          "gain_10p","loss_10p","gain_10q","loss_10q",
          "gain_11p","loss_11p","gain_11q","loss_11q",
          "gain_12p","loss_12p","gain_12q","loss_12q",
          "gain_13p","loss_13p","gain_13q","loss_13q",
          "gain_14p","loss_14p","gain_14q","loss_14q",
          "gain_15p","loss_15p","gain_15q","loss_15q",
          "gain_16p","loss_16p","gain_16q","loss_16q",
          "gain_17p","loss_17p","gain_17q","loss_17q",
          "gain_18p","loss_18p","gain_18q","loss_18q",
          "gain_19p","loss_19p","gain_19q","loss_19q",
          "gain_20p","loss_20p","gain_20q","loss_20q",
          "gain_21p","loss_21p","gain_21q","loss_21q",
          "gain_22p","loss_22p","gain_22q","loss_22q")
  
  file1=paste0("E://comp_allcancer//comp_",cancer,"//NDTinput//purity_",cancer,".txt")
  purity<-read.table(file1,header=T,sep = "\t")
  load("C://comp_allcancer//",cancer,"//timeall.rdata")
  
  setwd(paste0("C://comp_allcancer//",cancer,"//samplegroup"))
  i=1
  for(i in 1:length(scan)){
    gene1<-scan[i]
    m=sum(1,i)
    for (m in sum(1,i):length(scan)){
      gene2<-scan[m]
      sample1<-timeall[which(timeall$Event.Name==gene1),1]
      sample2<-timeall[which(timeall$Event.Name==gene2),1]
      sample<-intersect(sample1,sample2)
      if(length(sample)>30){
        sampleAB<-NA
        sampleBA<-NA
        sampleunkown<-NA
        j=1
        for(j in 1:length(sample)){
          tim<-timeall[timeall$Patient_ID==sample[j],]
          tim_1<-min(tim[which(tim$Event.Name==gene1),13])
          tim_2<-min(tim[which(tim$Event.Name==gene2),13])
          #Sometimes PhylogicNDT cannot infer an accurate false time probability density curve,result a mean of 0.05,removing these interference values
          if (tim_1!=0.05&tim_2!=0.05&tim_1<tim_2) {
            sample_name <-"sampleAB"
            sampleAB<-c(sampleAB,sample[j])
          }else if (tim_1!=0.05&tim_2!=0.05&tim_1>tim_2){
            sample_name <-"sampleBA"
            sampleBA<-c(sampleBA,sample[j])
          }else{
            sample_name <-"sampleAB_BA_unkown"
            sampleunkown<-c(sampleunkown,sample[j])
          }
        }
        sampleAB<-sampleAB[-1]
        sampleBA<-sampleBA[-1]
        sampleunkown<-sampleunkown[-1]
        
        sample1<-timeall[which(timeall$Event.Name==gene1),1]
        sample2<-timeall[which(timeall$Event.Name==gene2),1]
        sampleall<-unique(timeall$Patient_ID)
        sample0A<-setdiff(sample1,c(sampleAB,sampleBA))
        sample0B<-setdiff(sample2,c(sampleAB,sampleBA))
        sample00<-setdiff(sampleall,c(sample0A,sample0B,sampleAB,sampleBA))
        
        if(length(sampleAB)>15&length(sampleBA)>15){
          file7=paste0(gene1,"_",gene2,"_sample_A_B.rdata")
          save(sampleAB,sampleBA,sampleunkown,file=file7)
        }
      }
    }
  }
}  

#正常样本fkpm
get_gtex<-function(cancer_tissue){
  setwd(paste0("C://comp_allcancer//",cancer))
  gtex.exp=fread("C://comp_LUAD//gtex_RSEM_gene_fpkm.gz",header = T, sep = '\t',data.table = F)
  gtex.phe=read.table("C://comp_LUAD//GTEX_phenotype.tsv",header=T,sep="\t")
  gtex.pro=fread("C://comp_LUAD//gencode.v22.annotation.gene.probeMap",header = T, sep = '\t',data.table = F)
  gtex.pro=gtex.pro[,c(1,2)]
  gtex.exp=merge(gtex.exp,gtex.pro,by.x ="sample",by.y = "id" )
  colnames(gtex.phe)=c("Sample","body_site_detail (SMTSD)","primary_site","gender","patient","cohort")
  gtex.phe.s=filter(gtex.phe,primary_site==cancer_tissue)
  rownames(gtex.phe.s)=gtex.phe.s$Sample
  x1=intersect(rownames(gtex.phe.s),colnames(gtex.exp))
  gtex.s=gtex.exp[,c("sample",x1)]
  gtex.exp=merge(gtex.pro,gtex.s,by.x ="id",by.y ="sample")
  gtex.s1=distinct(gtex.exp,gene,.keep_all = T)
  gtex.s1 <- column_to_rownames(gtex.s1,"gene")
  gtex.s1 =gtex.s1[,-1]
  gtex.s1=2^gtex.s1
  gtex.s5=log2(gtex.s1-0.001+1)
  save(gtex.s5,file="gtex.s5_log2(fpkm+1).rdata")
}

#肿瘤样本fkpm
get_cancer_fpkm<-function(cancer){
  setwd(paste0("C://comp_allcancer//",cancer))
  XenaGenerate(subset = XenaHostNames=="gdcHub") %>% 
    XenaFilter(filterDatasets = "fpkm.tsv") %>% 
    XenaFilter(filterDatasets = cancer) -> df_todo
  XenaQuery(df_todo) %>%
    XenaDownload() -> xe_download
  expr = XenaPrepare(xe_download)
  expr.pro=data.table::fread("C://comp_LUAD//gencode.v22.annotation.gene.probeMap",header = T, sep = '\t',data.table = F)
  expr.pro=expr.pro[,c(1,2)]
  expr=merge(expr.pro,expr,by.y ="Ensembl_ID",by.x = "id" )
  expr=distinct(expr,gene,.keep_all = T)
  expr <- column_to_rownames(expr,"gene")
  expr<-expr[,which(substr(colnames(expr),14,15)=="01")]
  colnames(expr)<-substr(colnames(expr),1,12)
  save(expr,file=paste0(cancer,"_fkpm.rdata"))
}

#ssgsea
ssgsea_hallmark<-function(cancer){
  setwd(paste0("C://comp_allcancer//",cancer))
  load("C://comp_LUAD//hallmarkall.rdata")
  genesets <- split(msig_t2g_1$gene_symbol, msig_t2g_1$gs_name)
  load(paste0("C://comp_allcancer//",cancer,"_fpkm.rdata"))
  load("gtex.s5_log2(fpkm+1).rdata")
  allexpr<-merge(gtex.s5,expr,by= 0)
  allexpr <- column_to_rownames(expr,"Row.names")
  ssGSEA_hallmark <- gsva(expr = as.matrix(allexpr),  
                          gset.idx.list =genesets,
                          method = 'ssgsea',kcdf = 'Gaussian',abs.ranking = TRUE)
  ssgsea_hallmarkall<- t(scale(t(ssGSEA_hallmark)))
  ssgsea_hallmarkall[ssgsea_hallmarkall< -2] <- -2
  ssgsea_hallmarkall[ssgsea_hallmarkall>2] <- 2
  normalization <- function(x){return((x-min(x))/(max(x)-min(x)))}  
  ssgsea_hallmarkall <- normalization(ssgsea_hallmarkall)
  NM_ssgsea_hallmark<-ssgsea_hallmarkall[,1:ncol(gtex.s5)]
  NT_ssgsea_hallmark<-ssgsea_hallmarkall[,(ncol(gtex.s5)+1):ncol(allexpr)]
  save(NM_ssgsea_hallmark,NT_ssgsea_hallmark,ssgsea_hallmarkall,file = paste0(cancer,"ssgsea_hallmark.Rdata"))
}
ssgsea_kegg<-function(cancer){
  setwd(paste0("C://comp_allcancer//",cancer))
  kegg_geneset1 <- qusage::read.gmt("C://comp_LUAD//kegg_hsa.gmt") #返回的是list
  load("gtex.s5_log2(fpkm+1).rdata")
  load(paste0("C://comp_allcancer//",cancer,"_fpkm.rdata"))
  allexpr<-merge(gtex.s5,expr,by= 0)
  allexpr <- column_to_rownames(allexpr,"Row.names")
  ssGSEA_kegg <- gsva(expr = as.matrix(allexpr), 
                         gset.idx.list =kegg_geneset1,
                         method = 'ssgsea',kcdf = 'Gaussian',abs.ranking = TRUE)
  ssgsea_keggall<- t(scale(t(ssGSEA_kegg)))
  ssgsea_keggall[ssgsea_keggall< -2] <- -2
  ssgsea_keggall[ssgsea_keggall>2] <- 2
  normalization <- function(x){return((x-min(x))/(max(x)-min(x)))}  
  ssgsea_keggall <- normalization(ssgsea_keggall)
  NM_ssgsea_network<-ssgsea_networkall[,1:ncol(gtex.s5)]
  NT_ssgsea_network<-ssgsea_networkall[,(ncol(gtex.s5)+1):ncol(allexpr)]
  save(NM_ssgsea_network,NT_ssgsea_network,ssgsea_networkall,file = "ssgsea_network.Rdata")
}
ssgsea_K_network<-function(cancer){
  setwd(paste0("C://comp_allcancer//",cancer))
  load("C://comp_LUAD//network_genelist.rdata")
  load("gtex.s5_log2(fpkm+1).rdata")
  load(paste0("C://comp_allcancer//",cancer,"_fpkm.rdata"))
  allexpr<-merge(gtex.s5,expr,by= 0)
  allexpr <- column_to_rownames(allexpr,"Row.names")
  ssGSEA_network <- gsva(expr = as.matrix(allexpr), 
                         gset.idx.list =network,
                         method = 'ssgsea',kcdf = 'Gaussian',abs.ranking = TRUE)
  ssgsea_networkall<- t(scale(t(ssGSEA_network)))
  ssgsea_networkall[ssgsea_networkall< -2] <- -2
  ssgsea_networkall[ssgsea_networkall>2] <- 2
  normalization <- function(x){return((x-min(x))/(max(x)-min(x)))}  
  ssgsea_networkall <- normalization(ssgsea_networkall)
  NM_ssgsea_network<-ssgsea_networkall[,1:ncol(gtex.s5)]
  NT_ssgsea_network<-ssgsea_networkall[,(ncol(gtex.s5)+1):ncol(allexpr)]
  save(NM_ssgsea_network,NT_ssgsea_network,ssgsea_networkall,file = "ssgsea_network.Rdata")
}
ssgsea_GO<-function(cancer){
  setwd(paste0("C://comp_allcancer//",cancer))
  go_geneset1 <- qusage::read.gmt("C://comp_LUAD//c5.go.v2023.2.Hs.symbols.gmt") 
  load("gtex.s5_log2(fpkm+1).rdata")
  load(paste0(cancer,"_fkpm.rdata"))
  allexpr<-merge(gtex.s5,expr,by= 0)
  allexpr <- column_to_rownames(allexpr,"Row.names")
  load(paste0("C://comp_allcancer//",cancer,"_fpkm.rdata"))
  ssGSEA_go<- gsva(expr = as.matrix(allexpr), 
                         gset.idx.list =go_geneset1,
                         method = 'ssgsea',kcdf = 'Gaussian',abs.ranking = TRUE)
  ssgsea_goall<- t(scale(t(ssGSEA_go)))
  ssgsea_goall[ssgsea_goall< -2] <- -2
  ssgsea_goall[ssgsea_goall>2] <- 2
  normalization <- function(x){return((x-min(x))/(max(x)-min(x)))}  
  ssgsea_goall <- normalization(ssgsea_goall)
  NM_ssgsea_go<-ssgsea_goall[,1:ncol(gtex.s5)]
  NT_ssgsea_go<-ssgsea_goall[,(ncol(gtex.s5)+1):ncol(allexpr)]
  save(NM_ssgsea_go,NT_ssgsea_go,ssgsea_goall,file = paste0(cancer,"_ssgsea_go.Rdata"))
  
}

#get key mutation order pairs
get_keypairs<-function(cancer){
  load(paste0("C://comp_allcancer//",cancer,"_ssgsea_hallmark.Rdata"))
  colnames(ssgsea_hallmarkall)<-substr(colnames(ssgsea_hallmarkall),1,12)
  load(paste0("C://comp_allcancer//",cancer,"_ssgsea_kegg.Rdata"))
  colnames(ssgsea_networkall)<-substr(colnames(ssgsea_keggall),1,12)
  load(paste0("C://comp_allcancer//",cancer,"_ssgsea_network.Rdata"))
  colnames(ssgsea_networkall)<-substr(colnames(ssgsea_networkall),1,12)
  load(paste0("C://comp_allcancer//",cancer,"_ssgsea_go.Rdata"))
  colnames(ssgsea_goall)<-substr(colnames(ssgsea_goall),1,12)

  setwd(paste0("C://comp_allcancer//",cancer,"//allcomp"))
  fileall<-list.files(paste0("C://comp_allcancer//",cancer,"//samplegroup"))
  compall<-substr(fileall,1,nchar(fileall)-17)
  
 i=1
 for(i in 1:length(compall)){
    comp<-compall[i]
    load(paste0("C://comp_allcancer//",cancer,"//samplegroup//",comp,"_sample_A_B.rdata"))
    #运行functiondif
    sampleAB_ssgsea_hallmark <-ssgsea_hallmarkall[,intersect(colnames(ssgsea_hallmarkall),sampleAB)]
    sampleBA_ssgsea_hallmark <-ssgsea_hallmarkall[,intersect(colnames(ssgsea_hallmarkall),sampleBA)]
    sampleAB_ssgsea_kegg <-ssgsea_keggall[,intersect(colnames(ssgsea_keggall),sampleAB)]
    sampleBA_ssgsea_kegg <-ssgsea_keggall[,intersect(colnames(ssgsea_keggall),sampleBA)]
    sampleAB_ssgsea_network <-ssgsea_networkall[,intersect(colnames(ssgsea_networkall),sampleAB)]
    sampleBA_ssgsea_network <-ssgsea_networkall[,intersect(colnames(ssgsea_networkall),sampleBA)]
    sampleAB_ssgsea_go <-ssgsea_goall[,intersect(colnames(ssgsea_goall),sampleAB)]
    sampleBA_ssgsea_go <-ssgsea_goall[,intersect(colnames(ssgsea_goall),sampleBA)]
    H_comp<-NA
    n=1
    for (n in 1:nrow(sampleAB_ssgsea_hallmark)){
      p1<-wilcox.test(as.numeric(sampleAB_ssgsea_hallmark[n,]), as.numeric(sampleBA_ssgsea_hallmark[n,]), alternative = "two.sided", var.equal = FALSE,exact=FALSE)$p.value
      p2<-wilcox.test(as.numeric(sample00_ssgsea_hallmark[n,]), as.numeric(sampleAB_ssgsea_hallmark[n,]), alternative = "two.sided", var.equal = FALSE,exact=FALSE)$p.value
      p3<-wilcox.test(as.numeric(sample00_ssgsea_hallmark[n,]), as.numeric(sampleBA_ssgsea_hallmark[n,]), alternative = "two.sided", var.equal = FALSE,exact=FALSE)$p.value
      if(p1<0.05&(p2<0.05|p3<0.05)){
        H_data<-data.frame(comp=comp,pathway=rownames(sampleAB_ssgsea_hallmark)[n],patheay_num=n)
        H_data$logFC<-log2(mean(sampleBA_ssgsea_hallmark[n,]))-log2(mean(sampleAB_ssgsea_hallmark[n,]))
        H_data$P<-p1
        H_comp<-rbind(H_comp,H_data)
      }
    }
    H_comp<-as.data.frame(H_comp)
    if(nrow(H_comp)>1){
      setwd(paste0("C://comp_allcancer//",cancer,"//allcomp"))
      folder<-comp
      if (!dir.exists(folder)) {
        dir.create(folder)}
      setwd(folder)
      H_comp<-H_comp[-1,]
      H_comp$Padjust<-p.adjust(H_comp$P, method = "BH")
      save(H_comp,file="H_comp.rdata")
    }
    
    K_comp<-NA
    n=1
    for (n in 1:nrow(sampleAB_ssgsea_kegg)){
      p1<-wilcox.test(as.numeric(sampleAB_ssgsea_kegg[n,]), as.numeric(sampleBA_ssgsea_kegg[n,]), alternative = "two.sided", var.equal = FALSE,exact=FALSE)$p.value
      p2<-wilcox.test(as.numeric(sample00_ssgsea_kegg[n,]), as.numeric(sampleAB_ssgsea_kegg[n,]), alternative = "two.sided", var.equal = FALSE,exact=FALSE)$p.value
      p3<-wilcox.test(as.numeric(sample00_ssgsea_kegg[n,]), as.numeric(sampleBA_ssgsea_kegg[n,]), alternative = "two.sided", var.equal = FALSE,exact=FALSE)$p.value
      if(p1<0.05&(p2<0.05|p3<0.05)){
        K_data<-data.frame(comp=comp,pathway=rownames(sampleAB_ssgsea_kegg)[n],patheay_num=n)
        K_data$logFC<-log2(mean(sampleBA_ssgsea_kegg[n,]))-log2(mean(sampleAB_ssgsea_kegg[n,]))
        K_data$P<-p1
        K_comp<-rbind(K_comp,K_data)
      }
    }
    K_comp<-as.data.frame(K_comp)
    if(nrow(K_comp)>1){
      setwd(paste0("C://comp_allcancer//",cancer,"//allcomp"))
      folder<-comp
      if (!dir.exists(folder)) {
        dir.create(folder)}
      setwd(folder)
      K_comp<-K_comp[-1,]
      K_comp$Padjust<-p.adjust(K_comp$P, method = "BH")
      save(K_comp,file="K_comp.rdata")
    }
    
    network_comp<-NA
    m=1
    for (m in 1:nrow(sampleAB_ssgsea_network)){
      p1<-wilcox.test(as.numeric(sampleAB_ssgsea_network[m,]), as.numeric(sampleBA_ssgsea_network[m,]), alternative = "two.sided", var.equal = FALSE,exact=FALSE)$p.value
      p2<-wilcox.test(as.numeric(sample00_ssgsea_network[m,]), as.numeric(sampleAB_ssgsea_network[m,]), alternative = "two.sided", var.equal = FALSE,exact=FALSE)$p.value
      p3<-wilcox.test(as.numeric(sample00_ssgsea_network[m,]), as.numeric(sampleBA_ssgsea_network[m,]), alternative = "two.sided", var.equal = FALSE,exact=FALSE)$p.value
      if(p1<0.05&(p2<0.05|p3<0.05)){
        N_data<-data.frame(comp=comp,pathway=rownames(sampleAB_ssgsea_network)[m],pathway_num=m)
        N_data$logFC<-log2(mean(sampleBA_ssgsea_network[m,]))-log2(mean(sampleAB_ssgsea_network[m,]))
        N_data$P<-p1
        network_comp<-rbind(network_comp,N_data)
      }
    }
    network_comp<-as.data.frame(network_comp)
    if(nrow(network_comp)>1){
      setwd(paste0("C://comp_allcancer//",cancer,"//allcomp"))
      folder<-comp
      if (!dir.exists(folder)) {
        dir.create(folder)}
      setwd(folder)
      network_comp<-network_comp[-1,]
      network_comp$Padjust<-p.adjust(network_comp$P, method = "BH")
      save(network_comp,file="network_comp.rdata")
    }
    
    G_comp<-NA
    k=1
    for (k in 1:nrow(sampleAB_ssgsea_go)){
      p1<-wilcox.test(as.numeric(sampleAB_ssgsea_go[k,]), as.numeric(sampleBA_ssgsea_go[k,]), alternative = "two.sided", var.equal = FALSE,exact=FALSE)$p.value
      p2<-wilcox.test(as.numeric(sample00_ssgsea_go[k,]), as.numeric(sampleAB_ssgsea_go[k,]), alternative = "two.sided", var.equal = FALSE,exact=FALSE)$p.value
      p3<-wilcox.test(as.numeric(sample00_ssgsea_go[k,]), as.numeric(sampleBA_ssgsea_go[k,]), alternative = "two.sided", var.equal = FALSE,exact=FALSE)$p.value
      if(p1<0.05&(p2<0.05|p3<0.05)){
        G_data<-data.frame(comp=comp,pathway=rownames(sampleAB_ssgsea_go)[k],pathway_num=k)
        G_data$logFC<-log2(mean(sampleBA_ssgsea_go[k,]))-log2(mean(sampleAB_ssgsea_go[k,]))
        G_data$P<-p1
        G_comp<-rbind(G_comp,G_data)
      }
    }
    G_comp<-as.data.frame(G_comp)
    if(nrow(G_comp)>1){
      setwd(paste0("C://comp_allcancer//",cancer,"//allcomp"))
      folder<-comp
      if (!dir.exists(folder)) {
        dir.create(folder)}
      setwd(folder)
      G_comp<-G_comp[-1,]
      G_comp$Padjust<-p.adjust(G_comp$P, method = "BH")
      save(G_comp,file="G_comp.rdata")
    }

 }
}
