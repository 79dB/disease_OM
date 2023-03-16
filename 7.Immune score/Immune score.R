###计算免疫得分
library(estimate)
library(data.table)
library(ggplot2)
library(ggsignif)
library(ggpubr)

################################################################################
##免疫得分计算

##读入去批次后的表达谱
setwd("G:\\Lasso数据\\Lasso回归\\GPL570")
data <- as.data.frame(fread("RMA去批次结果.csv"))
info <- as.data.frame(fread("RMA样本信息.csv"))
##筛选四期患者
data <- data[,c(1,(1+which(colnames(data)[-1]%in%info$Sample)))]
##写出txt文件
setwd("G:\\Lasso数据\\免疫得分\\文件处理\\生出txt格式的表达谱")
#write.table(data,"txt格式表达谱.txt",sep="\t",col.names=T,row.names=F,quote=F)

##读入数据
data <- read.table("txt格式表达谱.txt",header=T)
##转换为gct格式
#filterCommonGenes(input.f="txt格式表达谱.txt", output.f="表达谱.txt", id="GeneSymbol")
##输出得分gct文件
#estimateScore(input.ds="表达谱.txt", output.ds="表达谱estimate_score.gct", platform="affymetrix")
##保存为其他格式
estimate_score <- read.table("表达谱estimate_score.gct", skip = 2, header = TRUE)

##写出csv格式文件
setwd("G:\\Lasso数据\\免疫得分\\免疫得分结果")
output <- estimate_score[,-2]
output <- t(output)
colnames(output) <- output[1,]
output <- output[-1,]
output <- as.data.frame(output)

##添加GSM编号和MRS分类
output$GSM_ID <- rownames(output)
##读入MRS信息
setwd("G:\\Lasso数据\\单因素和多因素")
info <- as.data.frame(fread("GEO数据的MRS信息.csv"))
##添加MRS
for (i in 1:nrow(output)) {
  output$MRS[i] <- info$Risk[info$GSM%in%output$GSM_ID[i]]
}
##调整顺序
output <- output[,c(5:6,1:4)]

setwd("G:\\Lasso数据\\免疫得分\\免疫得分结果")
##输出原格式结果
#write.csv(estimate_score,"原始免疫得分结果.csv",row.names = FALSE)
##输出转秩后列名为免疫得分的结果
#write.table(output,"免疫得分结果.csv",sep=",",col.names = T,row.names = F)

################################################################################
###绘制免疫得分相关箱型图
setwd("G:\\Lasso数据\\免疫得分\\免疫得分结果")
output <- as.data.frame(fread("免疫得分结果.csv"))
##基质的相关信息
emt <- output[,1:3]
##改变MRS
for (i in 1:nrow(emt)) {
  if(emt$MRS[i]%in%"MRS-high"){
    emt$MRS[i] <- "high-MRS"
  }else if(emt$MRS[i]%in%"MRS-low"){
    emt$MRS[i] <- "low-MRS"
  }
}
emt$MRS <- factor(emt$MRS,levels = c("high-MRS","low-MRS"))
#labels = c("high risk","low risk"))
##将得分转换为数值型
emt$StromalScore <- as.numeric(emt$StromalScore)
summary(emt$StromalScore)
##标准化
emt$StromalScore <- scale(emt$StromalScore)[,1]
summary(emt$StromalScore)

##绘制ggplot箱型图
p1 <- ggplot(data = emt,aes(x=MRS,y=StromalScore))+
  geom_boxplot(#outlier.shape = NA,##不显示离群点
    position=position_dodge(0.9),width=0.4,##调整箱图间距
    color=c("#CD5C5C","#1874CD"))+
  stat_boxplot(geom = "errorbar",width=0.3,color=c("#CD5C5C","#1874CD"))+  ##添加须线
  theme_bw()+
  theme(legend.position = "bottom",
        legend.key.size = unit(0.2, "inches"),
        panel.border = element_rect(),
        panel.grid = element_blank(),
        legend.text = element_text(size=12),
        axis.text = element_text(size=12,colour = "black"),
        axis.title = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        plot.title = element_text(size=12,hjust = 0.5),
        plot.margin = unit(c(0.2,0.8,0.2,0.1),"cm"))+
  xlab("")+ylab("Stromal score")+
  #scale_y_continuous(limits = c(-4,4))+
  ##添加显著性
  geom_signif(comparisons = list(c("high-MRS","low-MRS")),
              map_signif_level = T,textsize = 5)+
  scale_y_continuous(limits = c(-3,3))
p1
##添加空隙
p1 <- p1+theme(plot.margin = unit(c(0.8,0.2,1.5,0.2),"cm"))
#5*5
##免疫相关信息
imm <- output[,c(1,2,4)]
##改变MRS
for (i in 1:nrow(imm)) {
  if(imm$MRS[i]%in%"MRS-high"){
    imm$MRS[i] <- "high-MRS"
  }else if(imm$MRS[i]%in%"MRS-low"){
    imm$MRS[i] <- "low-MRS"
  }
}
imm$MRS <- factor(imm$MRS,levels = c("high-MRS","low-MRS"))
#labels = c("high risk","low risk"))
##将得分转换为数值型
imm$ImmuneScore <- as.numeric(imm$ImmuneScore)
summary(imm$ImmuneScore)
##标准化
imm$ImmuneScore <- scale(imm$ImmuneScore)[,1]
summary(imm$ImmuneScore)

##绘制ggplot箱型图
p2 <- ggplot(data = imm,aes(x=MRS,y=ImmuneScore))+
  geom_boxplot(#outlier.shape = NA,##不显示离群点
    position=position_dodge(0.9),width=0.4,##调整箱图间距
    color=c("#CD5C5C","#1874CD"))+
  stat_boxplot(geom = "errorbar",width=0.3,color=c("#CD5C5C","#1874CD"))+  ##添加须线
  theme_bw()+
  theme(legend.position = "bottom",
        legend.key.size = unit(0.2, "inches"),
        panel.border = element_rect(),
        panel.grid = element_blank(),
        legend.text = element_text(size=12),
        axis.text = element_text(size=12,colour = "black"),
        axis.title = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        plot.title = element_text(size=12,hjust = 0.5),
        plot.margin = unit(c(0.2,0.8,0.2,0.1),"cm"))+
  xlab("")+ylab("Immune score")+
  #scale_y_continuous(limits = c(-600,3000))+
  ##添加显著性
  geom_signif(comparisons = list(c("high-MRS","low-MRS")),
              annotations = "ns.",
              map_signif_level = T,textsize = 3)+
  scale_y_continuous(limits = c(-3,3))
p2
##添加空边
p2 <- p2+theme(plot.margin = unit(c(0.8,0.2,1.5,0.2),"cm"))
#5*5
##ESTIMATE得分
est <- output[,c(1,2,5)]
est$MRS <- factor(est$MRS,levels = c("MRS-high","MRS-low"))
#labels = c("high risk","low risk"))
##将得分转换为数值型
est$ESTIMATEScore <- as.numeric(est$ESTIMATEScore)
summary(est$ESTIMATEScore)
##标准化
est$ESTIMATEScore <- scale(est$ESTIMATEScore)[,1]
summary(est$ESTIMATEScore)

##绘制ggplot箱型图
p3 <- ggplot(data = est,aes(x=MRS,y=ESTIMATEScore))+
  geom_boxplot(#outlier.shape = NA,##不显示离群点
    position=position_dodge(0.9),width=0.4,##调整箱图间距
    color=c("#CD5C5C","#1874CD"))+
  stat_boxplot(geom = "errorbar",width=0.3,color=c("#CD5C5C","#1874CD"))+  ##添加须线
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10,colour = "black"),
        axis.title = element_text(size=14,colour = "black"),
        axis.text.x = element_text(hjust = 0.5,vjust = 1))+
  xlab("")+ylab("Estimate score")+
  #scale_y_continuous(limits = c(-600,3000))+
  ##添加显著性
  geom_signif(comparisons = list(c("MRS-high","MRS-low")),
              map_signif_level = T,textsize = 3)+
  theme(plot.margin = unit(c(1,1,0.5,1),"cm"))
p3

##肿瘤纯度
tum <- output[,c(1,2,6)]
tum$MRS <- factor(tum$MRS,levels = c("MRS-high","MRS-low"))
#labels = c("high risk","low risk"))
##将得分转换为数值型
tum$TumorPurity <- as.numeric(tum$TumorPurity)
summary(tum$TumorPurity)
##标准化
tum$TumorPurity <- scale(tum$TumorPurity)[,1]
summary(tum$TumorPurity)

##绘制ggplot箱型图
p4 <- ggplot(data = tum,aes(x=MRS,y=TumorPurity))+
  geom_boxplot(#outlier.shape = NA,##不显示离群点
    position=position_dodge(0.9),width=0.4,##调整箱图间距
    color=c("#CD5C5C","#1874CD"))+
  stat_boxplot(geom = "errorbar",width=0.3,color=c("#CD5C5C","#1874CD"))+  ##添加须线
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10,colour = "black"),
        axis.title = element_text(size=14,colour = "black"),
        axis.text.x = element_text(hjust = 0.5,vjust = 1))+
  xlab("")+ylab("Tumor purity")+
  #scale_y_continuous(limits = c(-600,3000))+
  ##添加显著性
  geom_signif(comparisons = list(c("MRS-high","MRS-low")),
              map_signif_level = T,textsize = 3)+
  theme(plot.margin = unit(c(1,1,0.5,1),"cm"))
p4


################################################################################
##免疫细胞类群分析

##读入网站分析结果
setwd("G:\\Lasso数据\\免疫得分\\免疫结果画箱型图")
data <- as.data.frame(fread("CIBERSORTx_Job1_Results.txt"))
##只保留L22相关信息
rownames(data) <- data$Mixture
data <- data[,-c(24:26)]
##更改列名
colnames(data) <- c("Mixture","B cell naive","B cell memory","B cell plasma",
                    "T cell CD8+","T cell CD4+ naive","T cell CD4+ memory resting",
                    "T cell CD4+ memory activated","T cell follicular helper",
                    "T cell regulatory (Tregs)","T cell gamma delta","NK cell resting",
                    "NK cell activated","Monocyte","Macrophage M0","Macrophage M1",
                    "Macrophage M2","Myeloid dendritic cell resting",
                    "Myeloid dendritic cell activated","Mast cell resting",
                    "Mast cell activated","Eosinophil","Neutrophil")
##列名
cname <- colnames(data)[-1]
rname <- rownames(data)
##读入样本MRS相关信息
setwd("G:\\Lasso数据\\单因素和多因素")
info <- as.data.frame(fread("GEO数据的MRS信息.csv")) 

##创建列表
l <- list()
for (i in 1:length(cname)) {
  l[[i]] <- data[,c(1,which(colnames(data)%in%cname[i]))]
}
##添加免疫细胞和MRS分类
for (i in 1:length(l)) {
  for (j in 1:nrow(l[[i]])) {
    l[[i]]$tissue <- colnames(l[[i]])[2]
  }
}
for (i in 1:length(l)) {
  for (j in 1:nrow(l[[i]])) {
    l[[i]]$MRS[j] <- info$Risk[info$GSM%in%l[[i]]$Mixture[j]]
  }
}
##调整顺序
for (i in 1:length(l)) {
  l[[i]] <- l[[i]][,c(1,4,3,2)]
}
##调整列名
for (i in 1:length(l)) {
  colnames(l[[i]])[4] <- "value"
}
##合并数据
result <- do.call(rbind,l)
##排序
result$tissue <- factor(result$tissue,levels = unique(result$tissue))
##查看顺序
levels(result$tissue)
##更改MRS
for (i in 1:nrow(result)) {
  if(result$MRS[i]%in%"MRS-high"){
    result$MRS[i] <- "high-MRS"
  }else if(result$MRS[i]%in%"MRS-low"){
    result$MRS[i] <- "low-MRS"
  }
}

##绘制箱型图
p <- ggplot(data = result,aes(x=tissue,y=value))+
  geom_boxplot(outlier.size = 0.5,aes(color=MRS),position=position_dodge(0.9),width=0.4)+
  stat_boxplot(mapping = aes(x=tissue,y=value,color=MRS),geom = "errorbar",width=0.3,position = position_dodge(0.9))+  ##添加须线
  scale_color_manual(values = c("#CD5C5C","#1874CD"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10,colour = "black"),
        axis.title = element_text(size=12,colour = "black"),
        axis.text.x = element_text(angle = 35,hjust = 1,vjust = 1))+
  xlab("")+ylab("Fraction")+
  #scale_y_continuous(limits = c(-4,4))+
  ##添加显著性
  geom_signif(comparisons = list(),
              map_signif_level = T)+
  theme(plot.margin = unit(c(1,1,0.5,1),"cm"))+
  stat_compare_means(aes(group=MRS),label = "p.signif",hide.ns = T,
                     label.y = c(0.15,0.2,0.28,0.28,0.13,0.33,0.2,0.2,0.2,0.2,0.18,
                                 0.13,0.15,0.3,0.3,0.32,0.25,0.2,0.15,0.4,0.15,0.4))
p$layers[[2]]$aes_params$textsize <- 5
p
#4/5*10
##统计NK激活的相关信息
tj <- result[result$tissue%in%"NK cell activated",]
summary(tj$value)
##查看不同MRS分类的中位数
low_med <- median(result$value[result$MRS%in%"MRS-low"])
low_med
up_med <- median(result$value[result$MRS%in%"MRS-high"])
up_med
