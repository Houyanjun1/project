library(ArrayExpress)
library(limma)
library(utils)
library(stringr)


library(ggVolcano)

vol_data <- read.table("Tumor-vs-Normal.limma.dif.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
#vol_data<-nrDEG[,c("logFC","P.Value")]
vol_data$row<-rownames(vol_data)
vol_adata <- add_regulate(vol_data, log2FC_name = "logFC",
                          fdr_name = "adj.P.Val",log2FC =1, fdr = 0.05)
table(vol_adata$regulate)
ggvolcano(vol_adata, x = "log2FoldChange", y = "padj",pointSize = 1,
          fills = c("#00a1d5","#efc000","#ee4c97"),
          colors = c("#00a1d5","#efc000","#ee4c97"),
          log2FC_cut = 1,
          FDR_cut = 0.05,
          add_line = TRUE,
          add_label = F,
          x_lab = "log2FoldChange",
          y_lab = "-log10pValue",
          label = "row",
          legend_position = "UL",
          output = FALSE)



setwd("F:\\20240904")  
library(pheatmap)          
inputFile="热图1.txt"      
groupFile="热图2.txt"  #分组文件   
outFile="heatmap.pdf"      
rt=read.table(inputFile,header=T,sep="\t",row.names=1,check.names=F)     #读取文件
ann=read.table(groupFile,header=T,sep="\t",row.names=1,check.names=F)    #读取样本属性文件


#绘制
pdf(file=outFile,width=7,height=5)
pheatmap(rt,
         annotation=ann,
         cluster_cols = F,
         cluster_rows = T,
         color = colorRampPalette(c("#00a1d5", "white", "#efc000"))(50),
         #color = colorRampPalette(c("blue", "white", "red"))(50),
         show_rownames=T,show_colnames=F,
         scale="row",  #矫正row
         #border_color ="NA",
         fontsize = 8,
         fontsize_row=7.5,
         fontsize_col=6)
dev.off()




#1LASSO
setwd("F:\\20240910")
library(glmnet)
train <- read.table("input.txt",sep="\t",header=T,check.names=F,row.names = 1)
train$group <- factor(train$group ,levels = c("0","1"))
dim(train)
train[1:4,1:4]

# 用LASSO-logitstic-Algorithm进行特征选择
# 转为lasso霃6?9要的格式
set.seed(20201015)
x <- as.matrix(train[,-1])  
(y <- ifelse(train$group == "0", 0,1)) #把分组信息换???01
#library(glmnet)


fit = glmnet(x, y, family = "binomial", nfold=10,alpha = 1, lambda = NULL)

# 画A???
pdf("lasso_CF.pdf", width = 6, height = 5)
plot(fit, xvar = "dev", label = TRUE)
dev.off()

cvfit = cv.glmnet(x, y, 
                  nfold=10, #例文描述???10-fold cross-validation
                  family = "binomial", type.measure = "class")
pdf("lasso_ME.pdf", width = 6, height = 5)
plot(cvfit)
dev.off()

cvfit$lambda.min #查看會6?9佳lambda

# 获取LASSO选出来的特征
myCoefs <- coef(cvfit, s="lambda.min");
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
(lasso_fea <- lasso_fea[-1])
# 把lasso找到的特征保存到文件
write.csv(lasso_fea,"feature_lasso.csv",row.names = F)

# 使用选出来的特征进行模型的预???
# 因为没有验证数据，这里用了训练数???
predict <- predict(cvfit, newx = x[1:nrow(x),], s = "lambda.min", type = "class")
table(predict,y)

#2SVM-RFE
library(e1071)
source("msvmRFE.R")
input <-read.table(file = "input.txt",header = T,sep = "\t",check.names=F,row.names = 1)
#采用十折交叉验证 (k-fold crossValidation???
# svmRFE(input, k = 10, halve.above = 50) #分割数据，分配随机数
svmRFE(input, k = 15, halve.above = 30) #分割数据，分配随机数
nfold = 15
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
# results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=50) #特征选择
results = lapply(folds, svmRFE.wrap, input, k=15, halve.above=30) #特征选择

top.features = WriteFeatures(results, input, save=F) #查看主要变量
head(top.features)
#把SVM-REF找到的特征保存到文件
write.csv(top.features,"feature_svm.csv",row.names = F)

# 运行时间主要取决于??6?9?择变量的个数，???6?9般的电脑还是不要选择太多变量
# 选前5个变量进行SVM模型构建，体验一???
featsweep = lapply(1:28, FeatSweep.wrap, results, input) #??һ???ܺ?ʱ???????÷?????????
featsweep


# 画图
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#dev.new(width=4, height=4, bg='white')
pdf("svm-error.pdf",width = 6,height = 5)
PlotErrors(errors, no.info=no.info) #查看错误???

dev.off()

# dev.new(width=4, height=4, bg='white')
pdf("svm-accuracy.pdf",width = 6,height = 5)
Plotaccuracy(1-errors,no.info=no.info) #查看准确???
dev.off()

which.min(errors) 
#3RF
library(tidyverse)
library(randomForest)
df <- read.table(file = "input.txt",
                 header = T,
                 sep = "\t",
                 row.names = 1)
head(df)[,1:5]
set.seed(1423)
df$group= as.factor(df$group)
fit.forest <- randomForest(group~., data=df,
                           ntree=3000)
fit.forest
pdf(file="forest.pdf",width = 6,height =5)
plot(fit.forest)
dev.off()
ch.min(fit.forest$err.rate[,1])

set.seed(1423)
fit.forest_2 <- randomForest(group~., data=df,
                             ntree=32)
fit.forest_2
pdf(file="forest_2.pdf",width = 6,height =6)
varImpPlot(fit.forest_2)
dev.off()
result_df <- data.frame(importance(fit.forest_2, type=2)) %>%
  arrange(desc(.))
write.table(result_df,file = "randomForest.txt",sep="\t")


setwd("F:\\20240508")
data3 <- read.table(file = "3.txt",header = T,sep = "\t",row.names = 1)#载入数据
library(regplot)
library(rms) 
dd <- datadist(data3)
options(datadist="dd")
f <- lrm(Group~ LYN+PLCG2+STAT5B+MMP9+IL6R, x=T, y=T, data=data3)#time.inc为时间增量
print(f)
nom <- nomogram(f, fun=plogis,
                lp=F, 
                fun.at = c(0.2,0.4,0.6,0.8,0.9,1.0),
                funlabel = "Risk"
)#执行Nomogram分析
par(mar = c(1, 3, 2, 2))
plot(nom,xfrac=.5,cex.axis=0.8,cex.var=1)
ckf
summary(f)

#建立校准曲线并绘制曲线图
cal1<-calibrate(f,method="boot",B=1000)
par(mar=c(6,5,1,2),cex = 1.0)#图形参数
plot(cal1,xlim=c(0,1.0),ylim=c(0,1.0),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),
     xlab = "Predicted Pr (pSS=1)", ylab = "Actual Probablity")



data  <- dca(data=data3, outcome="Group", predictors=c("predmodelA", "predmodelB","predmodelC"),smooth="TRUE", probability=c("TRUE", "TRUE","TRUE")) 

library(rmda)
data(data3)
set.seed(123)
BGN<- decision_curve(Group~ BGN,data = data3, 
                        family = binomial(link ='logit'),#模型类型，这里是二分类
                        thresholds= seq(0,1, by = 0.01),
                        confidence.intervals =0.95,#95可信区间
                        study.design = 'cohort')#研究类型，这里是队列研究
FN1<- decision_curve(Group~ FN1,data = data3, family = binomial(link ='logit'),
                        thresholds= seq(0,1, by = 0.01),
                        confidence.intervals =0.95,study.design ='cohort')
MMP2<- decision_curve(Group~ MMP2,data = data3, family = binomial(link ='logit'),
                      thresholds= seq(0,1, by = 0.01),
                      confidence.intervals =0.95,study.design ='cohort')
TP53<- decision_curve(Group~ TP53,data = data3, family = binomial(link ='logit'),
                      thresholds= seq(0,1, by = 0.01),
                      confidence.intervals =0.95,study.design ='cohort')
OSM<- decision_curve(Group~ OSM,data = data3, family = binomial(link ='logit'),
                      thresholds= seq(0,1, by = 0.01),
                      confidence.intervals =0.95,study.design ='cohort')
MMP9<- decision_curve(Group~ MMP9,data = data3, family = binomial(link ='logit'),
                       thresholds= seq(0,1, by = 0.01),
                       confidence.intervals =0.95,study.design ='cohort')
Nomogram<- decision_curve(Group  ~ LYN+PLCG2+STAT5B+MMP9+IL6R,data = data3,
                          family = binomial(link='logit'),
                          thresholds= seq(0,1, by = 0.01),
                          confidence.intervals =0.95,study.design ='cohort')

List<-list(BGN,FN1,MMP2,TP53,Nomogram)
plot_decision_curve(List,curve.names= c('BGN','FN1','MMP2','TP53','Nomogram'),
                    cost.benefit.axis =T,col = c('#00a1d5','#ee4c97','#efc000','#42b540','#925e9f'),
                    confidence.intervals =FALSE,standardize = F,
                    legend.position = "topright")#legend.position = "none"

plot_clinical_impact(Nomogram,
                     cost.benefit.axis = T,
                     n.cost.benefits= 8,
                     col =c('#ee4c97','#00a1d5'),
                     confidence.intervals= T,
                     ylim=c(0,1000),
                     legend.position="topright")





