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
specific<-PPI_table

PPI_Network <- graph.data.frame(specific,directed=FALSE)

PPI_Network <- igraph::simplify(PPI_Network, remove.multiple = TRUE, remove.loops = TRUE)

PPI_MultiplexObject <- create.multiplex(PPI_Network,Layers_Name=c("PPI"))

AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)

AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)

#########

seedgene<-read.table("DEG&driver.txt",sep="\t",header=F)
seedgene1<-as.matrix(seedgene)
##############

driver<-read.csv("cancer_gene_census.csv")
driver<-driver$Gene.Symbol
DEG2<-read.table("DEG2.txt",sep="\t",header=TRUE)
DEG2<-DEG2$Gene.Symbol
Gene<-unique(c(NPPI$SYMBOL.x,NPPI$SYMBOL.y))
seedgene<-setdiff(Gene, driver)
seedgene<-setdiff(seedgene,DEG2)

set.seed(123)
Res<-NULL
for (num in 1:100) {

seedgene1<-sample(seedgene,50)
#########


j=1
Result<-list()
while(j<=50){
  n<-which(specific[,1:2]==seedgene1[j])
  if(length(n)>0)
  {
    SeedGene<-seedgene1[j]
    RWR_PPI_Results<-Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,PPI_MultiplexObject,SeedGene)
    result<-RWR_PPI_Results[[1]][which(RWR_PPI_Results[[1]][2]>0.001),] 

    Result<-c(Result,list(c(SeedGene,t(result$NodeNames))))
  }
  j=j+1
}
###################

Disjoint<-0
Inter<-0
Contain<-0
for(i in 1:49){
  for (j in (i+1):50) {
    k<-length(intersect(Result[[i]],Result[[j]]))
    gene1<-length(Result[[i]])
    gene2<-length(Result[[j]])
    MIN<-min(gene1,gene2)
    if(k==0){
      Disjoint<-Disjoint+1
    }else if(k<MIN){
      Inter<-Inter+1
    }else{
      Contain<-Contain+1
    }
  }
}
Res<-rbind(Res,c(Disjoint,Inter,Contain))
}
setwd("G:/")
colnames(Res)<-c("Disjoint","Intersect","Contain")
write.table(Res,"Result.txt",quote=FALSE,sep="\t",row.names=FALSE)
#############

