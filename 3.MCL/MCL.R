#BiocManager::install("impute")
library(impute)
library(WGCNA)
library(ggplot2)
library(WeightedCluster)
F<-read.delim("G:/colon_cancer/Result/Sim_result_0.001.txt",header = FALSE)
G<-matrix(data=F$V3,nrow=46,ncol=46,byrow=T)

library(MCL)


TOM = TOMsimilarity(G)
dissTOM = 1-TOM

i<-1
Result<-NULL
for(i in seq(1,1.9, by=0.05)){
  
  res<-mcl(G,inflation =i,addLoops=T,max.iter = 1000)
  
  qual <- wcClusterQuality(dissTOM, res$Cluster)
  
  Result<-rbind(Result,c(i,qual$stats[5],qual$stats[4]))
}

colnames(Result)<-c("inflation","ASWw","ASW")
Result<-as.data.frame(Result)
ggplot(Result, aes(x=inflation, y=ASWw)) + 
  geom_line() + 
  geom_point(size=1)+
  annotate(geom = "point",
           x = 1.3, y = 0.22861986,color="red",size=3)+
  annotate(geom = "line",
           x =1.3, y = c(0,0.22),linetype="dashed")+
  scale_x_continuous(breaks=c(1, 1.1, 1.2, 1.3, 1.4, 1.5,1.6,1.7,1.8))+
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x = "inflation factor",y="Average Silhouette width (weighted)")
ggsave("G:/4.tiff", dpi=600)
dev.off()
