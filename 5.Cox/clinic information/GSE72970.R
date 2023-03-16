library(GEOquery)
library(data.table)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
setwd("F:\\R程序\\结直肠癌\\表达谱数据\\GEO数据\\RMA\\获取数据集样本临床信息\\GSE72970")
data <- as.data.frame(fread("GSE72970.csv",header = F))
##将数据框格式进行转换
sample <- as.character(data[2,-1])
data <- t(data)
colnames(data) <- data[1,]
data <- data[-1,]
rownames(data) <- data[,2]
data <- as.data.frame(data)
##保留相关信息
data2 <- data[,c(2,11:12,14:25)]
##患者相关信息
for (i in 1:ncol(data2)) {
  colnames(data2)[i] <- unique(unlist(strsplit(data2[,i],": "))[1])
}
colnames(data2)[1] <- "GSE_ID"
##将数据中的标题部分去掉
for (i in 1:nrow(data2)) {
  for (j in 2:(ncol(data2)-1)) {
    data2[i,j] <- unlist(strsplit(data2[i,j],": "))[2]
  }
}
##写出文件
setwd("G:\\Lasso数据\\GPL570临床信息整理新\\GSE72970")
write.table(data2,"GSE72970临床信息.csv",sep=",",col.names=T,row.names=F)
