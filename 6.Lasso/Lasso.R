library(data.table)
library(survival)
library(glmnet)
library(ggplot2)
library(grid)
library(timeROC)
library(RColorBrewer)
library(survivalROC)
library(survminer)
library(pheatmap)

##读入去批次后的表达谱
setwd("G:\\Lasso数据\\Lasso回归\\GPL570")
data <- as.data.frame(fread("RMA去批次结果.csv"))
info <- as.data.frame(fread("RMA样本信息.csv"))
##筛选四期患者
data <- data[,c(1,(1+which(colnames(data)[-1]%in%info$Sample)))]
##转秩
data <- t(data)
dim(data)
colnames(data) <- data[1,]
data <- data[-1,]
##转换数据框
data <- as.data.frame(data)

##读入88个基因
setwd("G:\\Lasso数据\\Lasso回归")
gene_info <- as.data.frame(fread("candidate gene.txt",header = F))
##查看基因交集,84个    
length(intersect(gene_info$V1,colnames(data)))
##查看未交上的基因
setdiff(gene_info$V1,colnames(data))

##保留交集基因的列
data <- data[,which(colnames(data)%in%intersect(gene_info$V1,colnames(data)))]

data$sample <- rownames(data)
##添加生存状态和生存时间
for (i in 1:nrow(data)) {
  data$status[i] <- info$overall_event[info$Sample%in%data$sample[i]]
  data$time[i] <- info$overall_survival[info$Sample%in%data$sample[i]]
}
data <- data[,c(85:87,1:84)]

##拆分生存时间的月份以及转换生存状态
for (i in 1:nrow(data)) {
  data$time[i] <- unlist(strsplit(data$time[i]," months"))
}
##将生存状态转换成数字型
data$time <- as.numeric(data$time)
#data$time <- as.integer(data$time*30)
for (i in 1:nrow(data)) {
  if(data$status[i]=="alive"){
    data$status[i] <- 0
  }else if(data$status[i]=="death"){
    data$status[i] <- 1
  }
}
data$status <- as.integer(data$status)
##将表达谱得分转化未数字型
for (i in 4:ncol(data)) {
  data[,i] <- as.numeric(data[,i])
}

##生存相关信息
survival_cancer <- data[,-1]
##获取基因
gene <- colnames(survival_cancer)[-c(1,2)]
##表达谱
x <- as.matrix(survival_cancer[,-c(1,2)])
##生存状态和时间
y <- survival_cancer[,c(2,1)]
names(y) <- c('time','status')
y$time <- as.double(y$time)
y$status <- as.double(y$status)
y <- as.matrix(Surv(y$time,y$status))

##LASSO回归图1
set.seed(5)
##下左上右
par(mai=c(0.8,0.8,0.5,0.3))
fit<-glmnet(x,y,family='cox')
plot(fit,label=F)
plot(fit,xvar="lambda",label=F,xlab="",cex.lab=2,cex.axis=2)
#,ylab="",xlab="")
#title(xlab= 'Log Lambda', ylab = 'Coefficients', line = 2)


##LASSO回归图2
par(mai=c(0.8,1,0.5,0.3))
lasso_fit <- cv.glmnet(x,y,family='cox',type.measure = 'deviance')
plot(lasso_fit)
#4*4/3*3

#拼接
# par(mai=c(1,1,0.6,0.3),oma=c(0.5,1,0.5,1))##下左上右
# par(mfrow=c(2,1))
# #par(mfcol=c(2,1),mar=c(1,4,1,1),oma=c(4,2,4,2))
# plot(lasso_fit,xlab="Log Lambda",xaxt = "n",yaxt = "n",cex.lab=2)#,xaxt="n")
# axis(1,-7:-2, cex.axis = 2,tcl=-0.3)
# axis(2,11:17, cex.axis = 2,tcl=-0.3)
# axis(3,c(81,80,81,80,78,73,66,53,41,29,22,6,1),cex.axis = 2)
# plot(fit,xvar="lambda",label=F,xaxt = "n",yaxt = "n",cex.lab=2)
# axis(1,-7:-2, cex.axis = 2,tcl=-0.3)
# axis(2,-3:3, cex.axis = 2,tcl=-0.3)
# #5*5

coefficient <- coef(lasso_fit,s=lasso_fit$lambda.min)
Active.Index <- which(as.numeric(coefficient)!=0)
active.coefficients <- as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox <- rownames(coefficient)[Active.Index]
##查看最低点处的值
length(sig_gene_multi_cox)
coef.min=as.matrix(coef(lasso_fit,s="lambda.min"))
#write.table(coef.min,"lasso_sig_gene(new).txt",sep="\t",col.names=FALSE,row.names=TRUE,quote=F)

##lasso系数绘图
dat_plot <- data.frame(Gene=sig_gene_multi_cox,Coefficient=active.coefficients)
dat_plot <- dat_plot[order(dat_plot$Coefficient),]

levels(dat_plot$Gene) <- dat_plot$Gene

dat_plot$Gene <- factor(dat_plot$Gene, levels=dat_plot$Gene, ordered=TRUE)
dat_plot$Group <- "high"
dat_plot[dat_plot$Coefficient<0,3] <- "low"
ggplot(data = dat_plot, mapping=aes(y = Coefficient, x = Gene, fill = Group)) + 
  ##stat="identity"表明取用样本点对应纵轴值
  # 条形图函数：position设置为"identity"是为了避免系统因绘制负值条形而引发的警告
  geom_bar(stat = 'identity', position = 'identity')+ 
  scale_fill_manual(values = c("#FFDDDD","#CCEEFF"))+
  ##横纵坐标互换
  coord_flip()+
  labs(y="LASSO coefficient")+
  ##transparent透明色
  theme(legend.position="none",panel.border=element_rect(fill='transparent'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 20,colour = "black"),
        axis.title.x = element_text(margin =margin(0.8,0,0,0,'cm')),
        axis.title.y = element_text(margin =margin(0,0.5,0,0,'cm')),
        axis.text = element_text(size = 20,colour = "black"),
        plot.margin = unit(c(1.5,1,0.5,1),"cm"))
##添加文字
grid.draw(grid.text(label="Protective",gp=gpar(col="blue",fontsize=20),x=unit(0.5,"npc"),y=unit(0.93,"npc")))
grid.draw(grid.text(label="Risk",gp=gpar(col="red",fontsize=20),x=unit(0.74,"npc"),y=unit(0.93,"npc")))
#5*5


##计算lasso风险得分
##表达谱
res_sur <- survival_cancer[,c(1:2,which(colnames(survival_cancer)%in%dat_plot$Gene))]
##系数
lasso_res <- dat_plot[,-3]
rownames(lasso_res) <- lasso_res$Gene
risk_score_table<-c()
for(i in 1:nrow(lasso_res)){
  gene<-rownames(lasso_res)[i]
  coef<-lasso_res$Coefficient[i]
  temp<-as.numeric(res_sur[,gene])*coef
  risk_score_table<-cbind(risk_score_table,temp)
}
colnames(risk_score_table)<-rownames(lasso_res)
rownames(risk_score_table)<-rownames(res_sur)
##添加风险得分，为每行的和
risk_score_table<-cbind(risk_score_table,'Risk_score'=rowSums(risk_score_table))
##为表达谱添加风险得分
risk_score_table<-cbind(res_sur,risk_score_table[,"Risk_score"])
colnames(risk_score_table)<-c(colnames(res_sur),"Risk_score")

##绘制ROC曲线
multi_ROC<-function(time_vector,risk_score_table){
  library('survivalROC')
  single_ROC<-function(single_time){
    ROC<-survivalROC(Stime=risk_score_table$time,
                     status=risk_score_table$status,
                     marker=risk_score_table$Risk_score,
                     predict.time=single_time,method='KM')
    data.frame('True_positive'=ROC$TP,'False_positive'=ROC$FP,
               'Cut_values'=ROC$cut.values,'Time_point'=rep(single_time,length(ROC$TP)),
               'AUC'=rep(ROC$AUC,length(ROC$TP)))
  }
  multi_ROC_list<-lapply(time_vector,single_ROC)
  do.call(rbind,multi_ROC_list)
}
#评估3-5年的AUC
for_multi_ROC<-multi_ROC(time_vector=c(12*seq(1,5,2)),risk_score_table)

##timeROC绘制ROC曲线
result<-with(risk_score_table,timeROC(T=time,
                                      delta=status,
                                      marker=Risk_score,
                                      cause=1,weighting="marginal",
                                      times=c(12*seq(1,5,2)),
                                      iid=TRUE))

dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(12,36,60)),each = nrow(result$TP)))


ggplot() + 
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,3,5),"-y: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=18,color="black"),
        axis.text = element_text(size=18,color="black"),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.65,0.18),
        legend.text = element_text(size=15),
        plot.margin =unit(c(1,1,1,1),"cm"),
        plot.title = element_text(size = 18,hjust = 0.5))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity",
       title = "Training (n = 257)")+
  coord_fixed()


#绘制AUC随时间变化
result <- with(risk_score_table,timeROC(T=time,
                                        delta=status,
                                        marker=Risk_score,
                                        cause=1,weighting="marginal",
                                        times=c(12*seq(1,5,0.1)),
                                        iid=TRUE))
par(mai=c(1,1,0.5,0.5))
plotAUCcurve(result, FP = 2, add = FALSE, conf.int = FALSE, conf.band = FALSE, col = "red")


##survivalROC确定最佳截点
##约登指数=灵敏度+特异度−1，即=TP−FP，约登指数最大时的截断值即为最佳截断值。
result <-with(risk_score_table, 
              survivalROC(Stime=time,
                          status=status,
                          marker=Risk_score,
                          predict.time=60,
                          method="KM"))
cutoff <- result$cut.values[which.max(result$TP-result$FP)]
cutoff
risk_score_table$Risk = ifelse(risk_score_table$Risk_score<cutoff,"MRS-low","MRS-high")
##计算生存曲线HR,log rank p
sdiff <- survdiff(Surv(time,status)~Risk,dat=risk_score_table)
p.val = 1 - pchisq(sdiff$chisq,length(sdiff$n)-1)
p.val
HR = (sdiff$obs[1]/sdiff$exp[1])/(sdiff$obs[2]/sdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdiff$exp[1]+1/sdiff$exp[2]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdiff$exp[1]+1/sdiff$exp[2])) ##计算km的HR
HR
low95
up95

##分组画KM曲线
fit <- survfit(Surv(time, status)~Risk, data=risk_score_table)
ggsurv<-ggsurvplot(fit, data=risk_score_table, pval=F,
                   #surv.median.line="hv", #添加中位生存时间线
                   palette = "jco",
                   risk.table = TRUE, 
                   risk.table.col="strata", #分线表颜色
                   risk.table.height=0.25,
                   ncensor.plot=FALSE,
                   NCENSOR.PLOT.HEIGHT=0.25,
                   legend.labs =c("high-MRS", "low-MRS"),  
                   conf.int = FALSE,xlab="Time in months",
                   break.x.by=12,xlim=c(0,60),
                   fontsize=8)
ggsurv$plot<- ggsurv$plot+labs(title="Training (OS)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 22),
        legend.title = element_text(size=22),
        legend.key.size=unit(0.5,'cm'))+
  annotate(geom="text", x=30, y=0.89, label="P = 4.404e-09\nHR = 2.556\n95%CI(1.918-3.405)",hjust=0,size=7)



#引入customize_labels 函数进行更高级的定制
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}

ggsurv <- customize_labels(
  ggsurv,
  font.title    = c(28, "bold", "darkblue"),#用长度为3的向量分别指定大小、类型、颜色
  font.subtitle = c(22, "bold", "purple"),
  font.caption  = c(22, "plain", "black"),
  font.x        = c(28, "bold", "darkred"),
  font.y        = c(28, "bold", "darkred"),
  font.xtickslab = c(24, "plain", "black"),
  font.ytickslab = c(24))
ggsurv

##热图
exp <- survival_cancer[,colnames(survival_cancer)%in%lasso_res$Gene]
##按照系数排序，将风险基因排在最前
exp <- exp[,as.character(lasso_res$Gene[order(lasso_res$Coefficient,decreasing = T)])]
risk_score_table <- risk_score_table[,c("status","time","Risk_score","Risk")]
new_data <- cbind(risk_score_table,exp)
table(new_data$Risk)
##按照risk_score升序排列
new_data_order <- new_data[order(new_data$Risk_score),]
##绘制热图
exp_order <- t(scale(new_data_order[,5:16]))
Risk <- data.frame(new_data_order[,"Risk"])
rownames(Risk) <- rownames(new_data_order)
colnames(Risk)="Risk"
colors=list(Risk=c("MRS-high"="red","MRS-low"="blue"))
par(mai=c(1,1,1,1))
#7*15
pheatmap(exp_order,
         show_clonames=F,
         scale="row",
         border=FALSE,
         frontsize_row=5,
         color=colorRampPalette(c("blue","white","red"))(256),
         show_rownames=T,
         show_colnames = F,
         cluster_rows=F,
         cluster_cols = F,
         legend_labels = "Expression",
         clustering_method = "average",
         cellheight=35,
         annotation_names_col = F,
         annotation_col = Risk,
         annotation_colors=colors)

##读入处理的TCGA二代数据(50个)
setwd("G:\\Lasso数据\\Lasso回归")
test <- as.data.frame(fread("TCGA测序四期患者.csv"))
rownames(test) <- test$sample
##将测试数据的生存时间转换为月
test$times <- test$times/30
colnames(test)[3] <- "time"

##进行检验
test <- test[,-1]
test <- test[,c(which(colnames(test)%in%colnames(res_sur)))]
##将表达谱转换为数字型
for (i in 3:ncol(test)) {
  test[,i] <- as.numeric(test[,i])
}
risk_score_table<-c()
for(i in 1:nrow(lasso_res)){
  gene<-rownames(lasso_res)[i]
  coef<-lasso_res$Coefficient[i]
  temp<-as.numeric(test[,gene])*coef
  risk_score_table<-cbind(risk_score_table,temp)
}
colnames(risk_score_table)<-rownames(lasso_res)
rownames(risk_score_table)<-rownames(test)
##添加风险得分，为每行的和
risk_score_table<-cbind(risk_score_table,'Risk_score'=rowSums(risk_score_table))
##为表达谱添加风险得分
risk_score_table<-cbind(test,risk_score_table[,"Risk_score"])
colnames(risk_score_table)<-c(colnames(test),"Risk_score")

##绘制ROC曲线
multi_ROC<-function(time_vector,risk_score_table){
  library('survivalROC')
  single_ROC<-function(single_time){
    ROC<-survivalROC(Stime=risk_score_table$time,
                     status=risk_score_table$status,
                     marker=risk_score_table$Risk_score,
                     predict.time=single_time,method='KM')
    data.frame('True_positive'=ROC$TP,'False_positive'=ROC$FP,
               'Cut_values'=ROC$cut.values,'Time_point'=rep(single_time,length(ROC$TP)),
               'AUC'=rep(ROC$AUC,length(ROC$TP)))
  }
  multi_ROC_list<-lapply(time_vector,single_ROC)
  do.call(rbind,multi_ROC_list)
}
#评估3-5年的AUC
for_multi_ROC<-multi_ROC(time_vector=c(12*seq(1,5,2)),risk_score_table)

##timeROC绘制ROC曲线
result<-with(risk_score_table,timeROC(T=time,
                                      delta=status,
                                      marker=Risk_score,
                                      cause=1,weighting="marginal",
                                      times=c(12*seq(1,5,2)),
                                      iid=TRUE))

dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(12,36,60)),each = nrow(result$TP)))


ggplot() + 
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,3,5),"-y: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=18,color="black"),
        axis.text = element_text(size=18,color="black"),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.65,0.18),
        legend.text = element_text(size=15),
        plot.margin =unit(c(1,1,1,1),"cm"),
        plot.title = element_text(size = 18,hjust = 0.5))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity",
       title = "Validation (n = 50)")+
  coord_fixed()


#绘制AUC随时间变化
result <- with(risk_score_table,timeROC(T=time,
                                        delta=status,
                                        marker=Risk_score,
                                        cause=1,weighting="marginal",
                                        times=c(12*seq(1,5,0.1)),
                                        iid=TRUE))
par(mai=c(1,1,0.5,0.5))
plotAUCcurve(result, FP = 2, add = FALSE, conf.int = FALSE, conf.band = FALSE, col = "red")

#######################################
##改变截点
##survivalROC确定最佳截点
##约登指数=灵敏度+特异度−1，即=TP−FP，约登指数最大时的截断值即为最佳截断值。
result <-with(risk_score_table, 
              survivalROC(Stime=time,
                          status=status,
                          marker=Risk_score,
                          predict.time=12,
                          method="KM"))
cutoff <- result$cut.values[which.max(result$TP-result$FP)]
cutoff
risk_score_table$Risk = ifelse(risk_score_table$Risk_score<cutoff,"MRS-low","MRS-high")
table(risk_score_table$Risk)

##计算生存曲线HR,log rank p
sdiff <- survdiff(Surv(time,status)~Risk,dat=risk_score_table)
p.val = 1 - pchisq(sdiff$chisq,length(sdiff$n)-1)
p.val
HR = (sdiff$obs[1]/sdiff$exp[1])/(sdiff$obs[2]/sdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdiff$exp[1]+1/sdiff$exp[2]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdiff$exp[1]+1/sdiff$exp[2])) ##计算km的HR
HR
low95
up95

##分组画KM曲线
fit <- survfit(Surv(time, status)~Risk, data=risk_score_table)
ggsurv<-ggsurvplot(fit, data=risk_score_table, pval=F,
                   #surv.median.line="hv", #添加中位生存时间线
                   palette = "jco",
                   risk.table = TRUE, 
                   risk.table.col="strata", #分线表颜色
                   risk.table.height=0.25,
                   ncensor.plot=FALSE,
                   NCENSOR.PLOT.HEIGHT=0.25,
                   legend.labs =c("high-MRS", "low-MRS"),  
                   conf.int = FALSE,xlab="Time in months",
                   break.x.by=12,xlim=c(0,60),
                   fontsize=8)
ggsurv$plot<-ggsurv$plot+labs(title="Validation (OS)")+theme(plot.title = element_text(hjust = 0.5),legend.text = element_text(size = 16),legend.title = element_text(size=16),legend.key.size=unit(0.5,'cm'))+
  annotate(geom="text", x=30, y=0.89, label="P = 1.385e-06\nHR = 6.135\n95%CI(1.713-21.972)",hjust=0,size=7)

#引入customize_labels 函数进行更高级的定制
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}

ggsurv <- customize_labels(
  ggsurv,
  font.title    = c(28, "bold", "darkblue"),#用长度为3的向量分别指定大小、类型、颜色
  font.subtitle = c(22, "bold", "purple"),
  font.caption  = c(22, "plain", "black"),
  font.x        = c(28, "bold", "darkred"),
  font.y        = c(28, "bold", "darkred"),
  font.xtickslab = c(24, "plain", "black"),
  font.ytickslab = c(24))
ggsurv

##热图
exp <- test[,colnames(test)%in%lasso_res$Gene]
##按照系数排序，将风险基因排在最前
exp <- exp[,as.character(lasso_res$Gene[order(lasso_res$Coefficient,decreasing = T)])]
risk_score_table <- risk_score_table[,c("status","time","Risk_score","Risk")]
new_data <- cbind(risk_score_table,exp)
table(new_data$Risk)
##按照risk_score升序排列
new_data_order <- new_data[order(new_data$Risk_score),]
##绘制热图
exp_order <- t(scale(new_data_order[,5:16]))
Risk <- data.frame(new_data_order[,"Risk"])
rownames(Risk) <- rownames(new_data_order)
colnames(Risk)="Risk"
colors=list(Risk=c("MRS-high"="red","MRS-low"="blue"))
##下左上右
par(mai=c(1,5,1,5))
##7*12
pheatmap(exp_order,
         show_clonames=F,
         scale="row",
         border=FALSE,
         frontsize_row=5,
         color=colorRampPalette(c("blue","white","red"))(256),
         show_rownames=T,
         show_colnames = F,
         cluster_rows=F,
         cluster_cols = F,
         legend_labels = "Expression",
         clustering_method = "average",
         cellheight=35,
         annotation_names_col = F,
         annotation_col = Risk,
         annotation_colors=colors)

