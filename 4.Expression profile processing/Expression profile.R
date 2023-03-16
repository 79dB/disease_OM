library(GEOquery)
library(data.table)

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

setwd("F:\\R程序\\结直肠癌\\表达谱数据\\GEO数据\\RMA\\获取数据集样本临床信息\\GSE17537")
gset <- getGEO('GSE17537', destdir=".",
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
clin <- clin[,c(2,10:18)]
##数据归位
clin[1:5,10] <- clin[1:5,9]
clin[1:5,9] <- clin[1:5,8]
clin[1:5,8] <- clin[1:5,7]
clin[1:5,7] <- clin[1:5,6]
clin[1:5,6] <- NA
##处理列名
for (i in 1:ncol(clin)) {
  colnames(clin)[i] <- unique(unlist(strsplit(clin[,i],": "))[1])
}
colnames(clin)[6] <- "grade"
colnames(clin) <- c("GSE_ID",colnames(clin)[-1])
##将数据中的标题部分去掉
for (i in 1:nrow(clin)) {
  for (j in 2:ncol(clin)) {
    clin[i,j] <- unlist(strsplit(clin[i,j],": "))[2]
  }
}
clin$tumor_location <- "colon"
##写出文件
setwd("G:\\Lasso数据\\GPL570临床信息整理新\\GSE17537")
write.table(clin,"GSE17537临床信息.csv",sep=",",col.names=T,row.names=F)

setwd("F:\\R程序\\结直肠癌\\表达谱数据\\GEO数据\\RMA\\获取数据集样本临床信息\\GSE29621")
gset <- getGEO('GSE29621', destdir=".",
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
clin <- pd[pd$characteristics_ch1.6%in%"ajcc staging: Stage 4",]
##保留相关信息
clin <- clin[,c(2,8,11:22)]
##处理列名
for (i in 1:ncol(clin)) {
  colnames(clin)[i] <- unique(unlist(strsplit(clin[,i],": "))[1])
}
colnames(clin)[1] <- "GSE_ID"
colnames(clin)[2] <- "tumor_location"
colnames(clin)[6] <- "m stage (0: no; 1: metastasis)"
##将数据中的标题部分去掉
for (i in 1:nrow(clin)) {
  for (j in 3:ncol(clin)) {
    clin[i,j] <- unlist(strsplit(clin[i,j],": "))[length(unlist(strsplit(clin[i,j],": ")))]
  }
}
##写出文件
setwd("G:\\Lasso数据\\GPL570临床信息整理新\\GSE29621")
write.table(clin,"GSE29621临床信息.csv",sep=",",col.names=T,row.names=F)

setwd("F:\\R程序\\结直肠癌\\表达谱数据\\GEO数据\\RMA\\获取数据集样本临床信息\\GSE39582")
gset <- getGEO('GSE39582', destdir=".",
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
clin <- pd[pd$characteristics_ch1.4%in%"tnm.stage: 4",]
##保留相关信息
clin <- clin[,c(2,8,12:42)]
##处理列名
for (i in 1:ncol(clin)) {
  colnames(clin)[i] <- unique(unlist(strsplit(clin[,i],": "))[1])
}
colnames(clin)[1] <- "GSE_ID"
colnames(clin)[2] <- "tumor_location"
clin$tumor_location <- "colorectal"
##将数据中的标题部分去掉
for (i in 1:nrow(clin)) {
  for (j in 3:ncol(clin)) {
    clin[i,j] <- unlist(strsplit(clin[i,j],": "))[length(unlist(strsplit(clin[i,j],": ")))]
  }
}
##去除os为0的样本
clin <- clin[-which(clin$`os.delay (months)`==0),]
##写出文件
setwd("G:\\Lasso数据\\GPL570临床信息整理新\\GSE39582")
write.table(clin,"GSE39582临床信息.csv",sep=",",col.names=T,row.names=F)

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

##样本信息
##读入临床信息
setwd("G:\\Lasso数据\\GPL570\\获取临床样本信息")
#GSE12945 <- as.data.frame(fread("GSE12945临床信息.csv"))
GSE17536 <- as.data.frame(fread("GSE17536临床信息.csv"))
GSE17537 <- as.data.frame(fread("GSE17537临床信息.csv"))
GSE29621 <- as.data.frame(fread("GSE29621临床信息.csv"))
GSE39582 <- as.data.frame(fread("GSE39582临床信息.csv"))
GSE72970 <- as.data.frame(fread("GSE72970临床信息.csv"))
#GSE83129 <- as.data.frame(fread("GSE83129临床信息.csv"))

#rownames(GSE12945) <- GSE12945$Symbol
rownames(GSE17536) <- GSE17536$Symbol
rownames(GSE17537) <- GSE17537$Symbol
rownames(GSE29621) <- GSE29621$Symbol
rownames(GSE39582) <- GSE39582$Symbol
rownames(GSE72970) <- GSE72970$Symbol
#rownames(GSE83129) <- GSE83129$Symbol
##查看样本交集
length(Reduce(intersect,list(#GSE12945$Symbol,
  GSE17536$Symbol,GSE17537$Symbol,
  GSE29621$Symbol,GSE39582$Symbol,GSE72970$Symbol)))
##删除样本
#GSE12945 <- GSE12945[,-ncol(GSE12945)]
GSE17536 <- GSE17536[,-ncol(GSE17536)]
GSE17537 <- GSE17537[,-ncol(GSE17537)]
GSE29621 <- GSE29621[,-ncol(GSE29621)]
GSE39582 <- GSE39582[,-ncol(GSE39582)]
GSE72970 <- GSE72970[,-ncol(GSE72970)]
#GSE83129 <- GSE83129[,-ncol(GSE83129)]
##添加平台
#GSE12945$GPL <- "GPL96"
GSE17536$GPL <- "GPL570"
GSE17537$GPL <- "GPL570"
GSE29621$GPL <- "GPL570"
GSE39582$GPL <- "GPL570"
GSE72970$GPL <- "GPL570"
#GSE83129$GPL <- "GPL6244"

##删除一些列
GSE39582 <- GSE39582[,c(1:8,10:13,31)]

##对数据进行操作，调整顺序
#GPL平台、年龄、性别、种族、癌症部位、转移数目、分期、辅助治疗\是否存活、复发、生存时间
#GSE12945 <- data.frame(Sample=rownames(GSE12945),
#                       GPL=GSE12945[,16],
#                       Age=GSE12945[,9],
#                       Sex=GSE12945[,8],
#                       ethnicity=NA,
#                       tumor_location=GSE12945[,1],
#                       synchronous_metastase=NA,
#                       meta_type=NA,
#                       number_of_meta=NA,
#                       pT=GSE12945[,11],
#                       pN=GSE12945[,12],
#                       pM=GSE12945[,13],
#                       Stage=GSE12945[,14],
#                       adjuvant_chemo=NA,
#                       overall_event=GSE12945[,4],
#                       dss_event=NA,
#                       dfs_event=NA,
#                       pfs_event=NA,
#                       rfs_event=NA,
#                       overall_survival=paste(GSE12945[,6],"months",sep=" "),
#                       dss_time=NA,
#                       dfs_time=NA,
#                       pfs_time=NA,
#                       rfs_time=NA)
#
GSE17536 <- data.frame(Sample=rownames(GSE17536),
                       GPL=GSE17536[,12],
                       Age=GSE17536[,1],
                       Sex=GSE17536[,2],
                       ethnicity=GSE17536[,3],
                       tumor_location=NA,
                       synchronous_metastase=NA,
                       meta_type=NA,
                       number_of_meta=NA,
                       pT=NA,
                       pN=NA,
                       pM=NA,
                       Stage=GSE17536[,4],
                       adjuvant_chemo=NA,
                       overall_event=GSE17536[,6],
                       dss_event=GSE17536[,7],
                       dfs_event=GSE17536[,8],
                       pfs_event=NA,
                       rfs_event=NA,
                       overall_survival=GSE17536[,9],
                       dss_time=GSE17536[,10],
                       dfs_time=GSE17536[,11],
                       pfs_time=NA,
                       rfs_time=NA)

GSE17537 <- data.frame(Sample=rownames(GSE17537),
                       GPL=GSE17537[,10],
                       Age=GSE17537[,1],
                       Sex=GSE17537[,2],
                       ethnicity=GSE17537[,3],
                       tumor_location=NA,
                       synchronous_metastase=NA,
                       meta_type=NA,
                       number_of_meta=NA,
                       pT=NA,
                       pN=NA,
                       pM=NA,
                       Stage=GSE17537[,4],
                       adjuvant_chemo=NA,
                       overall_event=GSE17537[,6],
                       dss_event=NA,
                       dfs_event=GSE17537[,7],
                       pfs_event=NA,
                       rfs_event=NA,
                       overall_survival=GSE17537[,8],
                       dss_time=NA,
                       dfs_time=GSE17537[,9],
                       pfs_time=NA,
                       rfs_time=NA)

GSE29621 <- data.frame(Sample=rownames(GSE29621),
                       GPL=GSE29621[,12],
                       Age=NA,
                       Sex=GSE29621[,1],
                       ethnicity=NA,
                       tumor_location=NA,
                       synchronous_metastase=NA,
                       meta_type=NA,
                       number_of_meta=NA,
                       pT=GSE29621[,2],
                       pN=GSE29621[,3],
                       pM=GSE29621[,4],
                       Stage=GSE29621[,6],
                       adjuvant_chemo=GSE29621[,7],
                       overall_event=GSE29621[,9],
                       dss_event=NA,
                       dfs_event=GSE29621[,11],
                       pfs_event=NA,
                       rfs_event=NA,
                       overall_survival=GSE29621[,8],
                       dss_time=NA,
                       dfs_time=GSE29621[,10],
                       pfs_time=NA,
                       rfs_time=NA)

GSE39582 <- data.frame(Sample=rownames(GSE39582),
                       GPL=GSE39582[,13],
                       Age=GSE39582[,2],
                       Sex=GSE39582[,1],
                       ethnicity=NA,
                       tumor_location=NA,
                       synchronous_metastase=NA,
                       meta_type=GSE39582[,7],
                       number_of_meta=NA,
                       pT=GSE39582[,4],
                       pN=GSE39582[,5],
                       pM=GSE39582[,6],
                       Stage=GSE39582[,3],
                       adjuvant_chemo=GSE39582[,8],
                       overall_event=GSE39582[,11],
                       dss_event=NA,
                       dfs_event=NA,
                       pfs_event=NA,
                       rfs_event=GSE39582[,9],
                       overall_survival=paste(GSE39582[,12],"months",sep=" "),
                       dss_time=NA,
                       dfs_time=NA,
                       pfs_time=NA,
                       rfs_time=GSE39582[,10])

GSE72970 <- data.frame(Sample=rownames(GSE72970),
                       GPL=GSE72970[,15],
                       Age=GSE72970[,2],
                       Sex=GSE72970[,1],
                       ethnicity=NA,
                       tumor_location=GSE72970[,4],
                       synchronous_metastase=GSE72970[,5],
                       meta_type=NA,
                       number_of_meta=NA,
                       pT=GSE72970[,7],
                       pN=GSE72970[,6],
                       pM=NA,
                       Stage=NA,
                       adjuvant_chemo=NA,
                       overall_event=GSE72970[,13],
                       dss_event=NA,
                       dfs_event=NA,
                       pfs_event=GSE72970[,11],
                       rfs_event=NA,
                       overall_survival=GSE72970[,14],
                       dss_time=NA,
                       dfs_time=NA,
                       pfs_time=GSE72970[,12],
                       rfs_time=NA)

##合并除GSE83129之外的结果
data <- do.call(rbind,list(#GSE12945,
  GSE17536,GSE17537,GSE29621,GSE39582,GSE72970))
##查看有无重复样本
length(data$Sample)
length(unique(data$Sample))

##统一生存信息
#os
for (i in 1:nrow(data)) {
  if(data$overall_event[i]%in%"1"){
    data$overall_event[i] <- "death"
  }else if(data$overall_event[i]%in%"0"){
    data$overall_event[i] <- "alive"
  }
}
for (i in 1:nrow(data)) {
  if(data$overall_event[i]%in%"dead"){
    data$overall_event[i] <- "death"
  }else if(data$overall_event[i]%in%"no death"){
    data$overall_event[i] <- "alive"
  }
}
#sex
for (i in 1:nrow(data)) {
  if(data$Sex[i]%in%"1"){
    data$Sex[i] <- "male"
  }else if(data$Sex[i]%in%"0"){
    data$Sex[i] <- "female"
  }
}
data$Sex <- tolower(data$Sex)
##添加转移个数
for (i in 1:nrow(data)) {
  if(data$pM[i]%in%c("1","M1")){
    data$number_of_meta[i] <- "1"
  }
}
##肿瘤部位归纳
for (i in 1:nrow(data)) {
  if(data$tumor_location[i]%in%c("caecum","Caecum")){
    data$tumor_location[i] <- "caecum"
  }else if(data$tumor_location[i]%in%c("colon flexure right","Left colon",
                                       "Right colon","Transverse colon","sigma")){
    data$tumor_location[i] <- "colon"
  }else if(data$tumor_location[i]%in%c("Rectum","Rectum-sigmoid junction",
                                       "rectum middle to lower third",
                                       "rectum upper third")){
    data$tumor_location[i] <- "rectum"
  }
}
for (i in 1:nrow(data)) {
  if(data$Stage[i]%in%c("4","IV","Stage 4")){
    data$Stage[i] <- "Stage 4"
  }
}
##去除盲肠的数据
#which(data$tumor_location=="caecum")
#data <- data[-which(data$tumor_location=="caecum"),]
##去除生存时间为0的数据
os <- data$overall_survival
os <- unlist(strsplit(os," months"))
os <- as.numeric(os)
which(os==0)
##GSE39582中的os=0
data <- data[-which(os==0),]

##写出文件
setwd("G:\\Lasso数据\\GPL570\\获取临床样本信息\\临床信息合并\\合并结果")
write.table(data,"RMA样本信息.csv",sep=",",col.names = T,row.names=F)

##表达谱去批次合并
##读入GPL570数据
setwd("G:\\Lasso数据\\GPL570")
GSE17536 <- as.data.frame(fread("GSE17536表达谱.csv"))
GSE17537 <- as.data.frame(fread("GSE17537表达谱.csv"))
GSE29621 <- as.data.frame(fread("GSE29621表达谱.csv"))
GSE39582 <- as.data.frame(fread("GSE39582表达谱.csv"))
GSE72970 <- as.data.frame(fread("GSE72970表达谱.csv"))
##列表
gpl570 <- list("GSE17536"=GSE17536,"GSE17537"=GSE17537,"GSE29621"=GSE29621,
               "GSE39582"=GSE39582,"GSE72970"=GSE72970)
##样本
#GPL570
sample_570 <- c()
for (i in 1:length(gpl570)) {
  sample_570 <- c(sample_570,colnames(gpl570[[i]])[-1])
}

##创建批次数据
#样本
sample_all <- c(sample_570)
#,sample_96)
length(sample_all)
length(unique(sample_all))
##创建数据框
info <- data.frame(sample_name=colnames(GSE17536)[-1],
                   batch=1)
info <- rbind(info,data.frame(sample_name=colnames(GSE17537)[-1],
                              batch=2))
info <- rbind(info,data.frame(sample_name=colnames(GSE29621)[-1],
                              batch=3))
info <- rbind(info,data.frame(sample_name=colnames(GSE39582)[-1],
                              batch=4))
info <- rbind(info,data.frame(sample_name=colnames(GSE72970)[-1],
                              batch=5))

info$sample_ID <- 1:nrow(info)
info <- info[,c(3,1,2)]

##读入拼接后的表达谱
setwd("G:\\Lasso数据\\GPL570\\表达谱合并\\合并结果")
data <- as.data.frame(fread("GEO数据RMA标准化合并结果.csv"))
rownames(data) <- data[,1]
data <- data[,-1]
# 提取batch的信息
batch = info$batch

# 去除批次效应
combat_edata1 = ComBat(dat=data, batch=info$batch, mod=NULL, par.prior=TRUE,  prior.plots=FALSE)
##转换数据框
combat_edata1 <- as.data.frame(combat_edata1)
combat_edata1$Symbol <- rownames(combat_edata1)
dim(combat_edata1)
combat_edata1 <- combat_edata1[,c(1007,1:1006)]

##写出文件
setwd("G:\\Lasso数据\\GPL570\\合并后去批次\\结果")
write.table(combat_edata1,"RMA去批次结果.csv",sep=",",col.names=T,row.names=F)

