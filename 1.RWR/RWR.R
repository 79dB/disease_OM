library(readxl)
#install.packages("BiocManager")
library("BiocManager")
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("RandomWalkRestartMH")
library(RandomWalkRestartMH)
library(igraph)
setwd("G:/colon_cancer")

PPI<-read_excel("41467_2019_9186_MOESM3_ESM.xlsx")
# transform id  
map1<- bitr(PPI$Protein_A_Entrez_ID, fromType = "ENTREZID",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
map2<- bitr(PPI$Protein_B_Entrez_ID, fromType = "ENTREZID",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
dt_merge <- merge(map1,PPI, by.y = "Protein_A_Entrez_ID", by.x = "ENTREZID")
PPI1<-merge(map2,dt_merge, by.y = "Protein_B_Entrez_ID", by.x = "ENTREZID")
PPI2<-PPI1[,c(2,4)]

NPPI<-PPI2[-which(PPI2$SYMBOL.x==PPI2$SYMBOL.y),]
####################################
PPI_table<-NPPI
seedgene<-read.table("FDEG&driver.txt",sep="\t",header=FALSE)
seedgene1<-as.matrix(seedgene)
specific<-PPI_table

PPI_Network <- graph.data.frame(specific,directed=FALSE)

PPI_Network <- igraph::simplify(PPI_Network, remove.multiple = TRUE, remove.loops = TRUE)

PPI_MultiplexObject <- create.multiplex(PPI_Network,Layers_Name=c("PPI"))

AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)

AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)

#########

setwd("G:/colon_cancer/FResult-0.001")
j=1
while(j<=46){
  n<-which(specific[,1:2]==seedgene1[j])
  if(length(n)>0)
  {
    SeedGene<-seedgene1[j]
    RWR_PPI_Results<-Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,PPI_MultiplexObject,SeedGene)
    result<-RWR_PPI_Results[[1]][which(RWR_PPI_Results[[1]][2]>0.001),] ##ãÐÖµ¿Éµ÷
    Result<-result


    assign(paste("result",j, sep=""),c(result$NodeNames,seedgene1[j]))
    
    assign(paste("Result",j, sep=""),c(Result$NodeNames))
    write.table(result,paste(seedgene1[j], ".txt", sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
    sep=t(c("Gene",SeedGene,dim(result)[1]))
    write.table(sep,"FRWR_result_0.001.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
  }
  j=j+1
}
###################

Disjoint<-NULL
Inter<-NULL
Contain<-NULL
for(i in 1:45){
  for (j in (i+1):46) {
    k<-length(intersect(get(paste("Result",i,sep = "")),get(paste("Result",j,sep = ""))))
    gene1<-length(get(paste("Result",i,sep = "")))
    gene2<-length(get(paste("Result",j,sep = "")))
    MIN<-min(gene1,gene2)
    if(k==0){
      Disjoint<-rbind(Disjoint,c(seedgene1[i],gene1,seedgene1[j],gene2,k))
    }else if(k<MIN){
      Inter<-rbind(Inter,c(seedgene1[i],gene1,seedgene1[j],gene2,k))
    }else{
      Contain<-rbind(Contain,c(seedgene1[i],gene1,seedgene1[j],gene2,k))
    }
  }
}
setwd("G:/colon_cancer/Result")
write.table(Disjoint,"FDisjoint.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(Inter,"FInter.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
#write.table(Contain,"FContain.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
#############

GENE<-c(t(NPPI$SYMBOL.x),t(NPPI$SYMBOL.y))
N<-length(unique(GENE))
i<-1
j<-2
for (i in 1:46) {
  a<-seedgene1[i]
  n<-length(get(paste("result",i,sep = "")))
  for (j in 1:46) {
    b<-seedgene1[j]
    m<-length(get(paste("result",j,sep = "")))
    k<-length(intersect(get(paste("result",i,sep = "")),get(paste("result",j,sep = ""))))
    p<-1-phyper(k-1,m,N-m,n)
    p.adjust<-p.adjust(p,method = "BH")
    sim<-length(intersect(get(paste("result",i,sep = "")),get(paste("result",j,sep = ""))))/
      length(union(get(paste("result",i,sep = "")),get(paste("result",j,sep = ""))))
    if(p.adjust<0.001){
      result=t(c(seedgene1[i],seedgene1[j],sim,p,p.adjust))
      write.table(result,"result_0.001.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
      sep=t(c(seedgene1[i],seedgene1[j],sim,p,p.adjust))
    }else{
      sep=t(c(seedgene1[i],seedgene1[j],0,p,p.adjust))
    }
    write.table(sep,"Sim_result_0.001.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
  }
}






