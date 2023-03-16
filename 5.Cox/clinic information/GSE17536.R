library(GEOquery)

setwd("F:\\R程序\\结直肠癌\\表达谱数据\\GEO数据\\RMA\\获取数据集样本临床信息\\GSE17536")
gset <- getGEO('GSE17536', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F)       ## 平台文件

class(gset)  #查看数据类型
length(gset)  #
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]] #
dat=exprs(a) #a现在是一个对象，取a这个对象通过看说明书知道要用exprs这个函数
dim(dat)#看一下dat这个矩阵的维度
dat[1:4,1:4] #查看dat这个矩阵的1至4行和1至4列，逗号前为行，逗号后为列
pd=pData(a) #通过查看说明书知道取对象a里的临床信息用pData

##筛选四期患者
clin <- pd[pd$characteristics_ch1.3%in%"ajcc_stage: 4",]
##保留相关信息
clin <- clin[,c(2,10:20)]
##处理列名
for (i in 1:ncol(clin)) {
  colnames(clin)[i] <- unique(unlist(strsplit(clin[,i],": "))[1])
}
colnames(clin) <- c("GSE_ID",colnames(clin)[-1])
##将数据中的标题部分去掉
for (i in 1:nrow(clin)) {
  for (j in 2:ncol(clin)) {
    clin[i,j] <- unlist(strsplit(clin[i,j],": "))[2]
  }
}
##添加肿瘤部位
clin$tumor_location <- "colon"

##写出文件
setwd("G:\\Lasso数据\\GPL570临床信息整理新\\GSE17536")
write.table(clin,"GSE17536临床信息.csv",sep=",",col.names=T,row.names=F)
