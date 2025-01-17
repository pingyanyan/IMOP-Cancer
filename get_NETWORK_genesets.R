#kegg子通路
library(KEGGgraph)
library(KEGGREST)
library(igraph)
kegg_geneset1 <- qusage::read.gmt("D://comp_LUAD//kegg_hsa.gmt") #返回的是list
pathway_all<-names(kegg_geneset1 )
pathway_id_all<-substr(pathway_all,1,8)
definition_allpathway<-list()
gene_allpathway<-list()

m=1
for(m in 1:length(pathway_all)){
  pathway_data <- keggGet(pathway_id_all[m])
  networkall<-pathway_data[[1]]$NETWORK$ELEMENT
  if(length(networkall)>0){
    geneall<-list()
    definitionall<-list()
    i=1
    for (i in 1:length(networkall)){
      network_data <- keggGet(substr(networkall[i],1,6))
      definition<-network_data[[1]][["DEFINITION"]][["DEFINITION"]]
      definitionall[[i]]<-definition
      gene_info<-network_data[[1]][["GENE"]]
      if(length(gene_info)>0){
        l<-length(gene_info)
        gene_ids<-gene_info[seq(1,l,2)]
        gene_names<-gene_info[seq(2,l,2)]
        gene_names<-unlist(lapply(gene_names,function(x) strsplit(x,split=";")[[1]][1]))
        geneall[[i]]<-gene_names
      }else{
        geneall[[i]]<-NA
      }
    }
    names(definitionall)<-networkall
    names(geneall)<-networkall
    definition_allpathway[[m]]<-definitionall
    gene_allpathway[[m]]<-geneall
  }else{
    definition_allpathway[[m]]<-NA
    gene_allpathway[[m]]<-NA
  }
}
names(definition_allpathway)<-pathway_id_all
names(gene_allpathway)<-pathway_id_all

save(definition_allpathway,gene_allpathway,file="network.rdata")

i=1
m=1
name<-NA
network<-list()
for (i in 1:length(gene_allpathway)){
  network1<-gene_allpathway[[i]]
  if(length(network1)>1||!is.na(network1)){
    j=1
    for (j in 1:length(network1)){
      network2<-gene_allpathway[[i]][[j]]
      network[[m]]<-network2
      name<-c(name,names(network1)[j])
      m=m+1
    }
  }
}
name<-name[-1]
names(network)<-name
save(network,file="network_genelist.rdata")
