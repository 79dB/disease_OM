library(data.table)
library(ggplot2)
library(grid)
library(timeROC)
library(RColorBrewer)
library(survivalROC)
library(survminer)
library(pheatmap)
library(forestplot)

##读入临床信息
setwd("G:\\Lasso数据\\GPL570临床信息整理新\\GSE17536")
GSE17536 <- as.data.frame(fread("GSE17536临床信息.csv"))
setwd("G:\\Lasso数据\\GPL570临床信息整理新\\GSE17537")
GSE17537 <- as.data.frame(fread("GSE17537临床信息.csv"))
setwd("G:\\Lasso数据\\GPL570临床信息整理新\\GSE29621")
GSE29621 <- as.data.frame(fread("GSE29621临床信息.csv"))
setwd("G:\\Lasso数据\\GPL570临床信息整理新\\GSE39582")
GSE39582 <- as.data.frame(fread("GSE39582临床信息.csv"))
setwd("G:\\Lasso数据\\GPL570临床信息整理新\\GSE72970")
GSE72970 <- as.data.frame(fread("GSE72970临床信息.csv"))
setwd("G:\\Lasso数据\\GPL570临床信息整理新\\TCGA二代临床数据新")
info <- as.data.frame(fread("TCGA临床数据.csv"))
##读入患者生存信息
setwd("F:\\R程序\\结直肠癌\\表达谱数据\\TCGA\\二代测序")
sur <- as.data.frame(fread("TCGA生存信息.txt"))

##读入MRS相关信息
setwd("G:\\Lasso数据\\单因素和多因素")
mrs_geo <- as.data.frame(fread("GEO数据的MRS信息.csv"))
mrs_tcga <- as.data.frame(fread("TCGA的MRS信息.csv"))

##调整数据框
GSE17536 <- data.frame(GSM_ID=GSE17536$GSE_ID,
                       GSE="GSE17536",
                       status=GSE17536$`overall_event (death from any cause)`,
                       os.time=GSE17536$`overall survival follow-up time`,
                       MRS=NA,
                       MRS_score=NA,
                       Age=GSE17536$age,
                       Gender=GSE17536$gender,
                       pT=NA,
                       pN=NA,
                       CEA_level=NA,
                       lymphatic_invasion=NA,
                       mutation=NA,   ##突变先写一行，有待确认具体数目
                       loss_expression_of_mismatch_repair_proteins=NA,
                       history_of_colon_polyps=NA,
                       preoperative_treatment=NA,
                       postoperative_treatment=NA,
                       meta_location=NA,
                       MSS=NA,
                       ly_tranfer=NA,
                       Peritoneal=NA)
##添加MRS
for (i in 1:nrow(GSE17536)) {
  GSE17536$MRS[i] <- mrs_geo$Risk[mrs_geo$GSM%in%GSE17536$GSM_ID[i]]
  GSE17536$MRS_score[i] <- mrs_geo$Risk_score[mrs_geo$GSM%in%GSE17536$GSM_ID[i]]
  
}

GSE17537 <- data.frame(GSM_ID=GSE17537$GSE_ID,
                       GSE="GSE17537",
                       status=GSE17537$`overall_event (death from any cause)`,
                       os.time=GSE17537$`overall survival follow-up time`,
                       MRS=NA,
                       MRS_score=NA,
                       Age=GSE17537$age,
                       Gender=GSE17537$gender,
                       pT=NA,
                       pN=NA,
                       CEA_level=NA,
                       lymphatic_invasion=NA,
                       mutation=NA,   ##突变先写一行，有待确认具体数目
                       loss_expression_of_mismatch_repair_proteins=NA,
                       history_of_colon_polyps=NA,
                       preoperative_treatment=NA,
                       postoperative_treatment=NA,
                       meta_location=NA,
                       MSS=NA,
                       ly_tranfer=NA,
                       Peritoneal=NA)
##添加MRS
for (i in 1:nrow(GSE17537)) {
  GSE17537$MRS[i] <- mrs_geo$Risk[mrs_geo$GSM%in%GSE17537$GSM_ID[i]]
  GSE17537$MRS_score[i] <- mrs_geo$Risk_score[mrs_geo$GSM%in%GSE17537$GSM_ID[i]]
}

GSE29621 <- data.frame(GSM_ID=GSE29621$GSE_ID,
                       GSE="GSE29621",
                       status=GSE29621$`os event`,
                       os.time=GSE29621$`overall survival (os)`,
                       MRS=NA,
                       MRS_score=NA,
                       Age=NA,
                       Gender=GSE29621$gender,
                       pT=GSE29621$`t stage`,
                       pN=GSE29621$`n stage`,
                       CEA_level=NA,
                       lymphatic_invasion=NA,
                       mutation=NA,   ##突变先写一行，有待确认具体数目
                       loss_expression_of_mismatch_repair_proteins=NA,
                       history_of_colon_polyps=NA,
                       preoperative_treatment=GSE29621$`preoperative chemo`,
                       postoperative_treatment=GSE29621$`adjuvant chemo`,
                       meta_location=NA,
                       MSS=NA,
                       ly_tranfer=NA,
                       Peritoneal=NA)
##添加MRS
for (i in 1:nrow(GSE29621)) {
  GSE29621$MRS[i] <- mrs_geo$Risk[mrs_geo$GSM%in%GSE29621$GSM_ID[i]]
  GSE29621$MRS_score[i] <- mrs_geo$Risk_score[mrs_geo$GSM%in%GSE29621$GSM_ID[i]]
}

##提取突变信息
mutation_geo582 <- GSE39582[,c("tp53.mutation","kras.mutation","braf.mutation")]

GSE39582 <- data.frame(GSM_ID=GSE39582$GSE_ID,
                       GSE="GSE39582",
                       status=GSE39582$os.event,
                       os.time=GSE39582$`os.delay (months)`,
                       MRS=NA,
                       MRS_score=NA,
                       Age=GSE39582$`age.at.diagnosis (year)`,
                       Gender=GSE39582$Sex,
                       pT=GSE39582$tnm.t,
                       pN=GSE39582$tnm.n,
                       CEA_level=NA,
                       lymphatic_invasion=NA,
                       mutation=NA,   ##突变先写一行，有待确认具体数目
                       loss_expression_of_mismatch_repair_proteins=GSE39582$mmr.status,
                       history_of_colon_polyps=NA,
                       preoperative_treatment=NA,
                       postoperative_treatment=GSE39582$chemotherapy.adjuvant,
                       meta_location=NA,
                       MSS=NA,
                       ly_tranfer=NA,
                       Peritoneal=NA)
##添加MRS
for (i in 1:nrow(GSE39582)) {
  GSE39582$MRS[i] <- mrs_geo$Risk[mrs_geo$GSM%in%GSE39582$GSM_ID[i]]
  GSE39582$MRS_score[i] <- mrs_geo$Risk_score[mrs_geo$GSM%in%GSE39582$GSM_ID[i]]
}


GSE72970 <- data.frame(GSM_ID=GSE72970$GSE_ID,
                       GSE="GSE72970",
                       status=GSE72970$`os censored`,
                       os.time=GSE72970$os,
                       MRS=NA,
                       MRS_score=NA,
                       Age=GSE72970$age,
                       Gender=GSE72970$Sex,
                       pT=GSE72970$pt,
                       pN=GSE72970$pn,
                       CEA_level=NA,
                       lymphatic_invasion=NA,
                       mutation=NA,   ##突变先写一行，有待确认具体数目
                       loss_expression_of_mismatch_repair_proteins=NA,
                       history_of_colon_polyps=NA,
                       preoperative_treatment="N",
                       postoperative_treatment=NA,
                       meta_location=NA,
                       MSS=NA,
                       ly_tranfer=NA,
                       Peritoneal=NA)
##添加MRS
for (i in 1:nrow(GSE72970)) {
  GSE72970$MRS[i] <- mrs_geo$Risk[mrs_geo$GSM%in%GSE72970$GSM_ID[i]]
  GSE72970$MRS_score[i] <- mrs_geo$Risk_score[mrs_geo$GSM%in%GSE72970$GSM_ID[i]]
}

##提取突变信息
mutation_tcga <- info[,c("kras_mutation_found")]
info$number_of_lymphnodes_positive_by_he <- as.numeric(info$number_of_lymphnodes_positive_by_he)
info$number_of_lymphnodes_positive_by_he[21] <- 10+info$number_of_lymphnodes_positive_by_he[21]
tcga <- data.frame(GSM_ID=info$sampleID,
                   GSE="TCGA",
                   status=NA,
                   os.time=NA,
                   MRS=NA,
                   MRS_score=NA,
                   Age=info$age_at_initial_pathologic_diagnosis,
                   Gender=info$gender,
                   pT=info$pathologic_T,
                   pN=info$pathologic_N,
                   CEA_level=info$preoperative_pretreatment_cea_level,
                   lymphatic_invasion=info$lymphatic_invasion,
                   kras.mutation=info$kras_mutation_found,   ##突变先写一行，有待确认具体数目
                   loss_expression_of_mismatch_repair_proteins=info$loss_expression_of_mismatch_repair_proteins_by_ihc,
                   history_of_colon_polyps=info$history_of_colon_polyps,
                   preoperative_treatment=info$history_of_neoadjuvant_treatment,
                   postoperative_treatment=info$postoperative_rx_tx,
                   meta_location=info$site_of_additional_surgery_new_tumor_event_mets,
                   MSS=info$CDE_ID_3226963,
                   ly_tranfer=info$number_of_lymphnodes_positive_by_he,
                   Peritoneal=info$pathologic_M)
##添加MRS
for (i in 1:nrow(tcga)) {
  tcga$MRS[i] <- mrs_tcga$Risk[mrs_tcga$GSM%in%tcga$GSM_ID[i]]
  tcga$MRS_score[i] <- mrs_tcga$Risk_score[mrs_tcga$GSM%in%tcga$GSM_ID[i]]
  tcga$status[i] <- sur$OS[sur$sample%in%tcga$GSM_ID[i]]
  tcga$os.time[i] <- sur$OS.time[sur$sample%in%tcga$GSM_ID[i]]/30
}


##为GEO数据添加突变信息
GSE17536 <- cbind(GSE17536[,1:12],data.frame(tp53.mutation=NA,kras.mutation=NA,
                                             braf.mutation=NA),
                  GSE17536[,14:21])
GSE17537 <- cbind(GSE17537[,1:12],data.frame(tp53.mutation=NA,kras.mutation=NA,
                                             braf.mutation=NA),
                  GSE17537[,14:21])
GSE29621 <- cbind(GSE29621[,1:12],data.frame(tp53.mutation=NA,kras.mutation=NA,
                                             braf.mutation=NA),
                  GSE29621[,14:21])
GSE39582 <- cbind(GSE39582[,1:12],mutation_geo582,GSE39582[,14:21])
GSE72970 <- cbind(GSE72970[,1:12],data.frame(tp53.mutation=NA,kras.mutation=NA,
                                             braf.mutation=NA),
                  GSE72970[,14:21])
tcga <- cbind(tcga[,1:12],data.frame(tp53.mutation=NA),tcga[,"kras.mutation"],
              data.frame(braf.mutation=NA),tcga[,14:21])
colnames(tcga) <- colnames(GSE17536)

##T、N分期转换为数字
for (i in 1:nrow(GSE39582)) {
  GSE39582$pT[i] <- unlist(strsplit(GSE39582$pT[i],"T"))[2]
  GSE39582$pN[i] <- unlist(strsplit(GSE39582$pN[i],"N"))[2]
}
##将加号转换为NA
GSE39582$pN[which(GSE39582$pN%in%"+")] <- NA

###
for (i in 1:nrow(GSE72970)) {
  GSE72970$pT[i] <- unlist(strsplit(GSE72970$pT[i],"pT"))[2]
  GSE72970$pN[i] <- unlist(strsplit(GSE72970$pN[i],"pN"))[2]
}
GSE72970$pT[which(GSE72970$pT%in%"X")] <- NA
GSE72970$pN[which(GSE72970$pN%in%"X")] <- NA

for (i in 1:nrow(tcga)) {
  tcga$pT[i] <- unlist(strsplit(tcga$pT[i],"T"))[2]
  tcga$pN[i] <- unlist(strsplit(tcga$pN[i],"N"))[2]
}
for (i in 1:nrow(tcga)) {
  if(tcga$pT[i]%in%c("4a","4b")){
    tcga$pT[i] <- 4
  }
  if(tcga$pN[i]%in%c("1a","1b","1c")){
    tcga$pN[i] <- 1
  }
  if(tcga$pN[i]%in%c("2a","2b")){
    tcga$pN[i] <- 2
  }
}

##拼接数据
data <- do.call(rbind,list(GSE17536,GSE17537,GSE29621,GSE39582,GSE72970))
##标准化MRS得分
data$MRS_score <- scale(data$MRS_score)
tcga$MRS_score <- scale(tcga$MRS_score)
data <- rbind(data,tcga)

##调整列名
colnames(data)[1] <- "Sample"
colnames(data)[2] <- "Source"
colnames(data)[6] <- "MRS_score"

##调整状态
for (i in 1:nrow(data)) {
  if(data$status[i]%in%c("0","alive","no death")){
    data$status[i] <- 0
  }else if(data$status[i]%in%c("1","dead","death")){
    data$status[i] <- 1
  }
}
##性别
data$Gender <- tolower(data$Gender)
##术前治疗
for (i in 1:nrow(data)) {
  if(data$preoperative_treatment[i]%in%c("N","No")){
    data$preoperative_treatment[i] <- "NO"
  }
}
##术后治疗
for (i in 1:nrow(data)) {
  if(data$postoperative_treatment[i]%in%c("Y","YES")){
    data$postoperative_treatment[i] <- "YES"
  }else if(data$postoperative_treatment[i]%in%c("N","NO")){
    data$postoperative_treatment[i] <- "NO"
  }else if(data$postoperative_treatment[i]%in%c("","N/A")){
    data$postoperative_treatment[i] <- NA
  }
}
##整理突变类型
for (i in 1:nrow(data)) {
  if(data$tp53.mutation[i]%in%"WT"){
    data$tp53.mutation[i] <- "NO"
  }else if(data$tp53.mutation[i]%in%"M"){
    data$tp53.mutation[i] <- "YES"
  }else if(data$tp53.mutation[i]%in%"N/A"){
    data$tp53.mutation[i] <- NA
  }
}
for (i in 1:nrow(data)) {
  if(data$kras.mutation[i]%in%"WT"){
    data$kras.mutation[i] <- "NO"
  }else if(data$kras.mutation[i]%in%"M"){
    data$kras.mutation[i] <- "YES"
  }else if(data$kras.mutation[i]%in%""){
    data$kras.mutation[i] <- NA
  }
}
for (i in 1:nrow(data)) {
  if(data$braf.mutation[i]%in%"WT"){
    data$braf.mutation[i] <- "NO"
  }else if(data$braf.mutation[i]%in%"M"){
    data$braf.mutation[i] <- "YES"
  }
}
##整理淋巴浸润
for (i in 1:nrow(data)) {
  if(data$lymphatic_invasion[i]%in%""){
    data$lymphatic_invasion[i] <- NA
  }
}
##纠错蛋白
for (i in 1:nrow(data)) {
  if(data$loss_expression_of_mismatch_repair_proteins[i]%in%""){
    data$loss_expression_of_mismatch_repair_proteins[i] <- NA
  }
}
##息肉史
for (i in 1:nrow(data)) {
  if(data$history_of_colon_polyps[i]%in%""){
    data$history_of_colon_polyps[i] <- NA
  }
}
##列名大写
cname <- colnames(data)
for (i in 1:length(cname)) {
  cname[i] <- paste(c(toupper(substr(cname[i],1,1)),
                      substr(cname[i],2,nchar(cname[i]))),collapse = "")
}
colnames(data) <- cname
##MMR类别修改
for (i in 1:nrow(data)) {
  if(data$Loss_expression_of_mismatch_repair_proteins[i]%in%"N/A"){
    data$Loss_expression_of_mismatch_repair_proteins[i] <- NA
  }else if(data$Loss_expression_of_mismatch_repair_proteins[i]%in%"dMMR"){
    data$Loss_expression_of_mismatch_repair_proteins[i] <- "YES"
  }else if(data$Loss_expression_of_mismatch_repair_proteins[i]%in%"pMMR"){
    data$Loss_expression_of_mismatch_repair_proteins[i] <- "NO"
  }
}
##是否肝转移
for (i in 1:nrow(data)) {
  if(data$Meta_location[i]%in%""){
    data$Meta_location[i] <- NA
  }else if(data$Meta_location[i]%in%"Liver"){
    data$Meta_location[i] <- "Yes"
  }else if(data$Meta_location[i]%in%"Other"){
    data$Meta_location[i] <- "No"
  }
}
##微卫星DNA,MSS和MSI
for (i in 1:nrow(data)) {
  if(data$MSS[i]%in%"MSI-H"|data$MSS[i]%in%"MSI-L"){
    data$MSS[i] <- "MSI"
  }
}
##淋巴转移
data$Ly_tranfer[which(data$Ly_tranfer>0)] <- "Yes"
data$Ly_tranfer[which(data$Ly_tranfer==0)] <- "No"

##腹膜转移
for (i in 1:nrow(data)) {
  if(data$Peritoneal[i]%in%"M1"){
    data$Peritoneal[i] <- "Unknow"
  }else if(data$Peritoneal[i]%in%c("M1a","M1b")){
    data$Peritoneal[i] <- "No"
  }
}


setwd("G:\\Lasso数据\\GPL570临床信息整理新\\GSE72970")
GSE72970 <- as.data.frame(fread("GSE72970临床信息.csv"))
##统计同时性转移
re1 <- GSE72970[,c(1,5)]
re2 <- info[,c(1,99)]
colnames(re2) <- colnames(re1)
##添加MRS分类
for (i in 1:nrow(re1)) {
  re1$MRS[i] <- mrs_geo$Risk[mrs_geo$GSM%in%re1$GSE_ID[i]]
}
for (i in 1:nrow(re2)) {
  re2$MRS[i] <- mrs_tcga$Risk[mrs_tcga$GSM%in%re2$GSE_ID[i]]
}


result <- rbind(re1,re2)
result <- result[,c(1,3,2)]
##统一格式
for (i in 1:nrow(result)) {
  if(result$`synchronous metastase`[i]%in%""){
    result$`synchronous metastase`[i] <- NA
  }else if(result$`synchronous metastase`[i]%in%"No"){
    result$`synchronous metastase`[i] <- "NO"
  }else if(result$`synchronous metastase`[i]%in%"Yes"){
    result$`synchronous metastase`[i] <- "YES"
  }
}
##输出带有同时性转移的数据
data$Synchronous <- NA
for (i in 1:nrow(data)) {
  if(length(result$`synchronous metastase`[result$GSE_ID%in%data$Sample[i]]!=0)){
    data$Synchronous[i] <- result$`synchronous metastase`[result$GSE_ID%in%data$Sample[i]]
  }
}

##查看分类
table(result$MRS,result$`synchronous metastase`)

##fisher检验
tj_syn <- as.matrix(table(result$MRS,result$`synchronous metastase`))
tj_syn
class(tj_syn)

fisher.test(tj_syn)

##写出文件
setwd("G:\\Lasso数据\\小修\\单因素和多因素森林图\\数据处理")
write.table(data,"Cox回归使用数据.csv",sep=",",col.names = T,row.names = F)



##绘制森林图
##读入临床表型信息
setwd("G:\\Lasso数据\\小修\\单因素和多因素森林图\\数据处理")
data <- as.data.frame(fread("Cox回归使用数据.csv"))
rownames(data) <- data$Sample
colnames(data)[3] <- "status"
colnames(data)[4] <- "time"

##所有数据进行单因素以及多因素生存分析
survival_cancer <- data[,-20]
##去除突变信息
survival_cancer <- survival_cancer[,-c(13:15)]
#13个表型作为单变量
##添加MRS，使用MRS值
variates <- colnames(survival_cancer)[c(6,7:20)]
variates
#方案2
uni_cox_in_bulk<-function(gene_list,survival_info_df){
  library('survival')
  uni_cox<-function(single_gene){
    formula<-as.formula(paste0('Surv(time,status)~',single_gene))
    surv_uni_cox<-summary(coxph(formula,data=survival_cancer))
    ph_hypothesis_p<-cox.zph(coxph(formula,data=survival_cancer))$table[1,3]
    {
      single_cox_report<-data.frame('uni_cox_sig_genes'=single_gene,
                                    'beta'=surv_uni_cox$coefficients[,1],
                                    'Hazard_Ratio'=surv_uni_cox$coefficients[,2],
                                    'HR.95L'=surv_uni_cox$conf.int[,3],
                                    'HR.95H'=surv_uni_cox$conf.int[,4],
                                    'pvalue'=surv_uni_cox$coefficients[,5])
      single_cox_report
    }
  }
  uni_cox_list<-lapply(variates,uni_cox)
  do.call(rbind,uni_cox_list)
}
uni_cox_df<-uni_cox_in_bulk(gene_list=variates,survival_info_df=survival_cancer)
##删除Inf
uni_cox_df <- uni_cox_df[uni_cox_df$HR.95H!="Inf",]
##处理名称
uni_cox_df$uni_cox_sig_genes <- c("MRS","Age","Gender","pT stage","pN stage","CEA level",
                                  "Lymphatic invasion","Loss expression of MMR",
                                  "History of colon polyps","Postoperative treatment",
                                  "Microsatellite status","Lymphnode metastasis",
                                  "Peritoneal metastasis","Metastatic presentation")
uni_cox_df <- uni_cox_df[c(1:5,11,14,6:7,12,8:10,13),]

#森林图可视化单因素cox分析结果
hz<-paste(round(uni_cox_df$Hazard_Ratio,3),
          "(",round(uni_cox_df$HR.95L,3),
          "-",round(uni_cox_df$HR.95H,3),")",sep="")
hz
tabletext<-cbind(c(NA,"Gene",uni_cox_df$uni_cox_sig_genes),
                 c(NA,"Coefficient",round(uni_cox_df$beta,3)),
                 c(NA,"P value",round(uni_cox_df$pvalue,3)),
                 c(NA,"Hazard Ratio(95% CI)",hz))
tabletext[35] <- "1.485e-14"
nrow(tabletext)

#5*8
forestplot(labeltext=tabletext,
           graph.pos=3,  #为Pvalue箱线图所在的位置
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           col=fpColors(box='#458B00', summary="#8B008B",lines="black", zero = "gray50"),
           mean=c(NA,NA,uni_cox_df$Hazard_Ratio),
           lower=c(NA,NA,uni_cox_df$HR.95L), #95%置信区间下限
           upper=c(NA,NA,uni_cox_df$HR.95H), #95%置信区间上限
           boxsize=0.4,#设置点估计的方形大小
           lwd.ci=2,   #设置区间估计线的粗细
           ci.vertices=TRUE, #置信区间的端点形状为T
           zero=1,#设置基准线x轴坐标
           lwd.zero=2, #设置基准线的粗细
           colgap=unit(5,"mm"),    #设置图形中的列间距
           xticks = c(0.5, 1,1.5), #横坐标刻度
           lwd.xaxis=2,            #设置X轴线的粗细
           lineheight = unit(1,"cm"), #设置图形中的行距
           graphwidth = unit(.3,"npc"), #图在表中的宽度比例
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
                           "17" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8),
                          xlab=gpar(cex = 1),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio",title="Univariate Cox regression analysis of clinical phenotype")

##对显著的表型接着进行多因素回归
survival_cancer <- data
survival_cancer <- survival_cancer[,-c(13:15)]
colnames(survival_cancer)[5] <- "Risk"
colnames(survival_cancer)[6] <- "MRS"

multi_cox<-coxph(Surv(time,status)~MRS+PT+PN,data=survival_cancer)
multi_cox<-summary(multi_cox)
gene<-c("MRS","PT","PN")
multi_cox_report<-data.frame('multi_cox_sig_genes'=gene,
                             'beta'=multi_cox$coefficients[,1],
                             'Hazard_Ratio'=multi_cox$conf.int[,1],
                             'HR.95L'=multi_cox$conf.int[,3],
                             'HR.95H'=multi_cox$conf.int[,4],
                             'pvalue'=multi_cox$coefficients[,5])
multi_cox_report$multi_cox_sig_genes <- c("MRS","pT stage","pN stage")
##更改p值
multi_cox_report$pvalue <- c("1.160e-10","0.137","0.113")


#森林图可视化多因素cox分析结果
hz<-paste(round(multi_cox_report$Hazard_Ratio,3),
          "(",round(multi_cox_report$HR.95L,3),
          "-",round(multi_cox_report$HR.95H,3),")",sep="")
tabletext<-cbind(c(NA,"Gene",multi_cox_report$multi_cox_sig_genes),
                 c(NA,"Coefficient",round(multi_cox_report$beta,3)),
                 c(NA,"P value",multi_cox_report$pvalue),
                 c(NA,"Hazard Ratio(95% CI)",hz))

nrow(tabletext)
#5*8
forestplot(labeltext=tabletext,
           graph.pos=3,  #为Pvalue箱线图所在的位置
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           col=fpColors(box='#458B00', summary="#8B008B",lines="black", zero = "gray50"),
           mean=c(NA,NA,multi_cox_report$Hazard_Ratio),
           lower=c(NA,NA,multi_cox_report$HR.95L), #95%置信区间下限
           upper=c(NA,NA,multi_cox_report$HR.95H), #95%置信区间上限
           boxsize=0.4,#设置点估计的方形大小
           lwd.ci=2,   #设置区间估计线的粗细
           ci.vertices=TRUE, #置信区间的端点形状为T
           zero=1,#设置基准线x轴坐标
           lwd.zero=2, #设置基准线的粗细
           colgap=unit(3,"mm"),    #设置图形中的列间距
           xticks = c(0.5, 1,1.5), #横坐标刻度
           lwd.xaxis=2,            #设置X轴线的粗细
           lineheight = unit(1,"cm"), #设置图形中的行距
           graphwidth = unit(.3,"npc"), #图在表中的宽度比例
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
                           "6" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex = 1),
                          title=gpar(cex = 1)),
           xlab="Hazard Ratio",title="Multivariate Cox regression")##training/validation


