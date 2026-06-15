library(ggpubr)
library(survival)
library(forestplot)
library(rms)
library(Hmisc)
library(rmda)
library(regplot)
library(ggplot2)
library(survcomp)
library(sva)
library(caret)
library(glmnet)
library(limma)
library(FactoMineR)
library(factoextra)
library(VennDiagram)
library(RColorBrewer)
library(ggrepel)
library(tidyverse)
library(pheatmap)

rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)      
rt[,"futime"]=rt[,"futime"]/365                            #Change the unit of survival time to years

#Loop through the groups and identify the groups that are significantly different for both the train and test sets.
for(i in 1:1000){
  #############Group the data#############
  inTrain<-createDataPartition(y=rt[,3],p=0.5,list=F)
  train<-rt[inTrain,]
  test<-rt[-inTrain,]
  trainOut=cbind(id=row.names(train),train)
  testOut=cbind(id=row.names(test),test)
  
  #############Univariate COX analysis#############
  outTab=data.frame()
  pFilter=0.05
  sigGenes=c("futime","fustat")
  for(i in colnames(train[,3:ncol(train)])){
    cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
    if(coxP<pFilter){
      sigGenes=c(sigGenes,i)
    }
  }
  train=train[,sigGenes]
  test=test[,sigGenes]
  uniSigExp=train[,sigGenes]
  uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
  
  #############Lasso Regression#############
  trainLasso=train
  trainLasso$futime[trainLasso$futime<=0]=0.003
  x=as.matrix(trainLasso[,c(3:ncol(trainLasso))])
  y=data.matrix(Surv(trainLasso$futime,trainLasso$fustat))
  fit <- glmnet(x, y, family = "cox", maxit = 1000)
  cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(coef != 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  lassoGene=c("futime","fustat",lassoGene)
  if(length(lassoGene)==2){
    next
  }	
  train=train[,lassoGene]
  test=test[,lassoGene]
  lassoSigExp=train
  lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
  
  #############Build a COX model#############
  multiCox <- coxph(Surv(futime, fustat) ~ ., data = train)
  multiCox=step(multiCox,direction = "both")
  multiCoxSum=summary(multiCox)
  
  #Output model-related information
  outMultiTab=data.frame()
  outMultiTab=cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
  
  #Output the risk file of the train group
  riskScore=predict(multiCox,type="risk",newdata=train)          
  coxGene=rownames(multiCoxSum$coefficients)
  coxGene=gsub("`","",coxGene)
  outCol=c("futime","fustat",coxGene)
  medianTrainRisk=median(riskScore)
  risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
  trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))
  
  #Output the risk files of the test group
  riskScoreTest=predict(multiCox,type="risk",newdata=test)      
  riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
  testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest)); 
  
  diff=survdiff(Surv(futime, fustat) ~risk,data = train)
  pValue=1-pchisq(diff$chisq,df=1)
  roc = survivalROC(Stime=train$futime, status=train$fustat, marker = riskScore, predict.time =1,  method="KM")
  
  diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
  pValueTest=1-pchisq(diffTest$chisq,df=1)
  rocTest = survivalROC(Stime=test$futime, status=test$fustat, marker = riskScoreTest, predict.time =1,  method="KM")
  
  if((pValue<0.05) & (roc$AUC>0.75) & (pValueTest<0.05) & (rocTest$AUC>0.75)){
    #Output grouped results
    write.table(trainOut,file="04.train.txt",sep="\t",quote=F,row.names=F)
    write.table(testOut,file="04.test.txt",sep="\t",quote=F,row.names=F)
    #Output the results of a single factor
    write.table(outTab,file="05.uniCox.xls",sep="\t",row.names=F,quote=F)
    write.table(uniSigExp,file="05.uniSigExp.txt",sep="\t",row.names=F,quote=F)
    #Output the results of Lasso regression
    write.table(lassoSigExp,file="06.lassoSigExp.txt",sep="\t",row.names=F,quote=F)
    pdf("06.lambda.pdf")
    plot(fit, xvar = "lambda", label = TRUE)
    dev.off()
    pdf("06.cvfit.pdf")
    plot(cvfit)
    abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
    dev.off()
    #Output multi-factor results
    write.table(outMultiTab,file="07.multiCox.xls",sep="\t",row.names=F,quote=F)
    write.table(testRiskOut,file="riskTest.txt",sep="\t",quote=F,row.names=F)
    write.table(trainRiskOut,file="riskTrain.txt",sep="\t",quote=F,row.names=F)
    break
  }
}

#Draw a forest map
options(forestplot_new_page = FALSE)
pdf(file="07.forest.pdf",width = 8,height = 5)
ggforest(multiCox,main = "Hazard ratio",cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()


###survival analysis####

inputFile="riskTest.txt"        
survFile="test-survival.pdf"        
rocFile="test-ROC.pdf"               

inputFile="riskTrain.txt"        
survFile="train-survival.pdf"         
rocFile="train-ROC.pdf"  


rt=read.table(inputFile,header=T,sep="\t")

#?Ƚϸߵͷ????????????죬?õ???????pֵ
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
#survival plot####
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=T,
                   pval=pValue,
                   pval.size=5,
                   risk.table=TRUE,
                   legend.labs=c("High risk", "Low risk"),
                   legend.title="Risk",
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c('#EBBA37',"#8AD293"),
                   risk.table.height=.25,
                   title = "Merged Test")
pdf(file=survFile,onefile = FALSE,width = 6.5,height =5.5)
print(surPlot)
dev.off()


###ROC####
ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
               marker=rt$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file=rocFile,width=7,height=5)
plot(ROC_rt,time=1,col="#8AD293",title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='#EBBA37',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col="#4665D9",add=TRUE,title=FALSE,lwd=2)
title(main = "Merged Test", cex.main = 1.5)
legend('bottomright',
       c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],3)),
         paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],3)),
         paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],3))),
       col=c("#8AD293",'#EBBA37',"#4665D9"),lwd=2,bty = 'n')
dev.off()

####Comparison of clinical information####
inputFile="TARGET_clinical.txt"#input file
outFile="TARGET-boxplot.pdf"      #output file
inputFile="GEO_clinical.txt"
outFile="GEO-boxplot.pdf"

#reading a file
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
x=colnames(rt)[2]
y=colnames(rt)[3]
colnames(rt)=c("id","Group","Expression")

#Set up a comparison group
group=levels(factor(rt$Group))
rt$Group=factor(rt$Group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#Draw a boxplot
boxplot=ggboxplot(rt, x="Group", y="Expression", color="Group",
                  xlab=x,
                  ylab=y,
                  legend.title=x,
                  palette = c("#F1515E","#1DBDE6"),
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons,method =  "wilcox.test")+labs(x='',y='Futime(years)',title = "TARGET")
#stat_compare_means(comparisons = my_comparisons,method =  "wilcox.test")+labs(x='',y='Futime(years)',title = "GES21257")
#save
pdf(file=outFile,width=6,height=5)
print(boxplot)
dev.off()

####multicox analysis####
rm(list = ls())
clrs <- fpColors(box="black",line="black", summary="black")            
rt=read.table("clinical.txt",header=T,sep="\t",check.names=F,row.names=1)

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)


rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          
pdf(file="forest-multicox.pdf",onefile = FALSE,
    width = 6,             
    height = 4,           
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()

####Nomogram and calibration curve####
#Read the risk input file
rt=read.table("clinical.txt", header=T, sep="\t", check.names=F, row.names=1)
dd <- datadist(rt)
options(datadist = "dd")

colnames(rt)

#model building

fit <- coxph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, data = rt)              # BpM

pdf(file="BpM_Nomogram.pdf", width=6, height=6)
nom=regplot(fit,
            clickable=F,
            title="",
            points=TRUE,
            droplines=T,
            observation=NULL,
            rank=NULL,
            failtime = c(1,3,5),
            showP = F,
            prfail = F) 
dev.off()
#calibration curves#
#1 year
pdf(file="BpMcalibration.pdf", width=6, height=6)
fit1 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(fit1, cmethod="KM", method="boot", u=1, m=28,B=1000)##  m = 28 indicates that 20% of the total sample size.
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="Firebrick2", sub=F)
#3 year
fit1 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(fit1, cmethod="KM", method="boot", u=3, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="MediumSeaGreen", sub=F, add=T)
#5 year
fit1 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(fit1, cmethod="KM", method="boot", u=5, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("Firebrick3","MediumSeaGreen","NavyBlue"), lwd=3, bty = 'n')
dev.off()


##remove batch####
rm(list = ls())
rt=read.table("mergegroup.txt",sep="\t",header=T,check.names=F)
rt <- t(rt)
rt=as.matrix(rt)

exp <- rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

expr <- data
batch <- c(rep("GSE21257",53),rep("TARGET",86))
tissue <- c(rep("Non_metastasis",19),rep("Metastasis",34),rep("Non_metastasis",65),rep("Metastasis",21))
mode <- model.matrix(~as.factor(tissue))

####remove Batch Effect
limma_expr <- removeBatchEffect(expr,batch = batch,design = mode)

#PCA analysis without removing batch effects
pre.pca <- PCA(t(expr),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= "point",
             title = "Before Batchremoving",
             col.ind = batch,
             addEllipses = TRUE,
             legend.title="Group"  )

#Batch Effect-Corrected PCA Analysis
combat.pca <- PCA(t(limma_expr),graph = FALSE)
fviz_pca_ind(combat.pca,
             geom= "point",
             title = "After Batchremoving",
             col.ind = batch,
             addEllipses = TRUE,
             legend.title="Group"  )

####venn plot####
plot <- venn.diagram(
  x = list(set1,set2,set3),
  category.names = c("ARG", "TARGET","GSE21257"),
  filename =NULL,  
  output=FALSE,
  # Circle Properties:
  col = "black", 
  lty = 1, 
  lwd = 1, 
  fill = c("#C1C1C1","#7AAACB","#E89874"),
  alpha = 0.60, 
  label.col = "black",
  cex = .5, 
  fontfamily = "serif",
  fontface = "bold",
  
  # Collection Name Attribute:
  cat.col = c("#C1C1C1","#7AAACB","#E89874"),
  cat.cex = .6,
  cat.fontfamily = "serif"
)

####Volcano plot####
#input file
rm(list = ls())
df <- read.csv("limmaTab.csv",header = T,row.names = 1)
head(df)

df$Group <- factor(ifelse(df$P.Value < 0.05 & abs(df$logFC) >= 0.2,
                          ifelse(df$logFC >= 0.2, 'Up','Down'),'Stable'))
df[1:10,1:7]

table(df$Group)
df$gene <- row.names(df)

p <- ggplot(df, aes(x = logFC, y = -log10(P.Value),,colour = Group))+
  geom_point( shape = 19, size=2.5,stroke = 0.5)+
  
  scale_color_manual(values=c( "#1874CD",'gray',"#CD2626"))+
  ylab('-log10 (Pvalue)')+
  xlab('log2 (Fold Change)')+
  labs(title = "No_metastases vs Metastases")+
  #The gene name of the added focus point
  geom_text_repel(
    data = df[df$P.Value < 0.05 & abs(df$logFC) > 0.2,],
    aes(label = gene),
    size = 3.5,
    segment.color = NA )+ 
  geom_vline(xintercept = c(-0.2,0.2),lty = 2, col = "black", lwd = 0.5)+
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black", lwd = 0.5)+
  theme_bw(
    base_line_size = 1  
  )+
  guides(fill = guide_legend(override.aes = list(size =3)))+
  theme_bw()+
  theme(
    axis.title.x = element_text(hjust = 0.5),
    legend.position = c(0.08, 0.86)
  )


p + theme(  panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(
              hjust = 0.5, 
              size = 14,   
              face = "bold", 
              vjust = 1.5
            )
)

####Pheatmap####

# Load the expression data and annotation data
DEG <- read.table("diffExp.txt",header=T, sep="\t", check.names=F,row.names = 1)

annotation_col <- read.table("annotation_col.txt",header=T, sep="\t", check.names=F,row.names = 1)

annotation_row <- read.table("annotation_row.txt",header=T, sep="\t", check.names=F,row.names = 1)
#Convert to factor type
annotation_col[] <- lapply(annotation_col, as.factor)

annotation_row[] <- lapply(annotation_row, as.factor)
#Set the name of the row for the comment
rownames(annotation_row) <- rownames(DEG)

rownames(annotation_col) <- colnames(DEG)

#Set annotation color
annotation_colors <- list(
  Regulation=c("Up"="#FC9F5B","Down"="#25998F"),
  Gender = c("Male"="#c9bc9c",
             "Female"="#4665d9"
  ),
  Fustat = c("Alive"="#ff91c2",
             "Dead"="#3e3a39"),
  Group = c("Metastases"="#911fb4","No_metastases"="#818001")
)

# Draw a heat map
p <- pheatmap(DEG,
              scale = "row",
              cluster_rows = F,
              cluster_cols = F,
              annotation_row = annotation_row,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              color = colorRampPalette(rev(brewer.pal(10, "RdBu")))(50),
              show_rownames = TRUE,
              show_colnames =FALSE,
              treeheight_row = 0,
              treeheight_col = 0 )
p

dev.off()
