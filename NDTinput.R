#Download the input file for integrating TCGA data into PhylogicNDT
#Download maf, ascat
get_maf<-function(cancer){
  query <- GDCquery(project = paste0("TCGA-",cancer), 
                    data.category = "Simple Nucleotide Variation",
                    data.type = "Masked Somatic Mutation",
                    access = "open")
  GDCdownload(query)
  file<-paste0("TCGA-",cancer,"_SNP.Rdata")
  GDCprepare(query, save = T,save.filename = file)
  query <- GDCquery(
  project = paste0("TCGA-",cancer),
  data.category = "Copy Number Variation",
  data.type = "Allele-specific Copy Number Segment",
  workflow.type = "ASCAT2")
  GDCdownload(query)
  filename<-paste0("TCGA-",cancer,"_ASCAT.Rdata")
  GDCprepare(query, save = T,save.filename = filename)
  ##read data
  a<-data[,c("Tumor_Sample_Barcode","Hugo_Symbol","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","HGVSp","t_ref_count","t_alt_count")]
  load(filename)
  b<-data[,c(2,3,4,6,7,8)]
  a$Tumor_Sample_Barcode<-substr(a$Tumor_Sample_Barcode, 1,12)
  a$Chromosome<-substring(a$Chromosome, 4)
  a$Hugo_Symbol<-paste(a$Hugo_Symbol,"_",a$HGVSp)
  allsample<-unique(a[,1])
  b <-b[which(substr(b$Sample,16,16) == "A"),]
  b$Chromosome<-substring(b$Chromosome, 4)
  b$Sample<-substr(b$Sample, 1,12)
  filename1<-paste0(cancer,"_maf.rdata")
  save(a,file=filename1)
  filename2<-paste0(cancer,"_ascat.rdata")
  save(b,file=filename2)
}

#Download clinical information
get_clinical<-function(cancer){
  query <- GDCquery(project = paste0("TCGA-",cancer),
                    data.category = "Clinical",
                    data.type = "Clinical Supplement",
                    data.format = "BCR XML")
  GDCdownload(query)
  clinical <- GDCprepare_clinic(query, clinical.info = "patient")
  save(clinical,file="clinical.Rdata")
}

#NDTinput
sampleinput<-function(cancer){
  file1<-paste0(cancer,"_maf.rdata")
  load(file1)
  file2<-paste0(cancer,"_ascat.rdata")
  load(file2)
  allsample<-unique(a[,1])
  allsample1<-unique(b[,6])
  allsample<-as.character(t(allsample))
  allsample1<-as.character(t(allsample1))
  allsample2<-intersect(allsample,allsample1)
  allsample2<-as.data.frame(allsample2)
  ##Match maf and ascat
  i=1
  for(i in 1:nrow(allsample2)){
    c<-subset(a, Tumor_Sample_Barcode == as.character(allsample2[i,]))
    d<-subset(b, Sample == as.character(allsample2[i,]))
    c<-c[,c(2,3,4,5,6,8,9)]
    d<-d[,c(1,2,3,4,5)]
    c$Chromosome = as.character(c$Chromosome)
    c$Start_Position = as.numeric(c$Start_Position)
    d$Chromosome = as.character(d$Chromosome)
    d$Start = as.numeric(d$Start)
    d$End = as.numeric(d$End)
    npos <- length(c$Start_Position)
    local_cn_a1<- vector(mode="numeric",length=npos)
    local_cn_a2<- vector(mode="numeric",length=npos)
    Chrs <- unique(c$Chromosome)
    for(chr in Chrs){
      t1.chr.matches <- which(c$Chromosome==chr)
      t2.chr.compares <- which(d$Chromosome==chr)
      for(match in t1.chr.matches){
        cn_a1=0
        cn_a2=0
        candidate = as.numeric(as.vector(c$Start_Position)[match])
        for(compare in t2.chr.compares){
          comp.start <- as.numeric(as.vector(d$Start)[compare])
          comp.stop <- as.numeric(as.vector(d$End)[compare])
          if(candidate>=comp.start&&candidate<=comp.stop){
            cn_a1 = as.numeric(as.vector(d$Minor_Copy_Number)[compare])
            cn_a2 = as.numeric(as.vector(d$Major_Copy_Number)[compare])
            break
          }
        }
        local_cn_a1[match] = cn_a1
        local_cn_a2[match] = cn_a2
      }
    }
    Tab1 <- as.data.frame(cbind(c,local_cn_a1,local_cn_a2))
    colnames(Tab1)<-c("Hugo_Symbol","Chromosome","Start_position","Reference_Allele","Tumor_Seq_Allele2","t_ref_count","t_alt_count","local_cn_a1","local_cn_a2")
    outfile<-paste(allsample2[i,],"_input",".txt",sep="")
    write.table(Tab1,file=outfile,quote = F,sep="\t",row.names = F)
  }
  ##purity
  purity<-read.table("D://PhylogicNDT//purity.txt",header=T,sep = "\t")
  rownames(purity)<-purity[,1]
  allsample2<-intersect(allsample,allsample1)
  purity1<-purity[allsample2,]
  purity1<-purity1[,c(1,4)]
  purity1<-na.omit(purity1)
  purity1$array<-substr(purity1$array,1,12)
  filename3=paste0("purity_",cancer,".txt")
  write.table(purity1,file=filename3,quote = F,sep="\t",row.names = F )
  allsample3<-rownames(purity1)
  save(allsample3,file="sample.rdata")
  ## integrate ASCAT input file 
  colnames(b)<-c("Chromosome","Start_position","End_Position","A2_CN","A1_CN","sample")
  b$Chromosome = as.character(b$Chromosome)
  b$Start_position = as.numeric(b$Start_position)
  b$End_Position = as.numeric(b$End_Position)
  b$A1_CN = as.numeric(b$A1_CN)
  b$A2_CN = as.numeric(b$A2_CN)
  b$sample = as.character(b$sample)
  
  allsample2<-as.data.frame(allsample2)
  i=1
  for(i in 1:nrow(allsample2)){
    d<-subset(b, sample == as.character(allsample2[i,]))
    d<-d[,c(1,2,3,5,4)]
    outfile<-paste(allsample2[i,],"_ascat",".txt",sep="")
    write.table(d,file=outfile,quote = F,sep="\t",row.names = F)
  }
}
