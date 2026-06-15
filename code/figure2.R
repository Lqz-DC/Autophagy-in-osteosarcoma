library(survival)
library(forestplot)
library(survivalROC)
library(rms)
library(Hmisc)
library(rmda)
library(regplot)
library(ggplot2)
library(survcomp)

####unicox####
clrs <- fpColors(box="green",line="darkblue", summary="royalblue")             
rt=read.table("input.txt",header=T,sep="\t",check.names=F,row.names=1)

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)


rt=read.table("uniCox.xls",header=T,sep="\t",row.names=1,check.names=F)
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
pdf(file="forest-uniCox.pdf",onefile = FALSE,
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


####multicox####
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")            
rt=read.table("input.txt",header=T,sep="\t",check.names=F,row.names=1)

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
pdf(file="forest-multi.pdf",onefile = FALSE,
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
####multiROC####
rt=read.table("input.txt",header=T,sep="\t",check.names=F,row.names=1)    
rt$futime=rt$futime/365
rocCol=rainbow(ncol(rt)-2)
aucText=c()


pdf(file="multiROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("riskScore"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)


j=1
for(i in colnames(rt[,3:(ncol(rt)-1)])){
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =1, method="KM")
  j=j+1
  aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()


####Nomogram and calibration diagram####
#Read the risk input file
rt=read.table("all.txt", header=T, sep="\t", check.names=F, row.names=1)
dd <- datadist(rt)
options(datadist = "dd")

#merging data
colnames(rt)
rt1=rt[,-c(6,7,8,9,11)]#Metastasis,Age,Gender and risk Score
rt2=rt1[,-6]#Metastasis,Age and Gender
rt3 <- rt1[,-5]#Age,Gender and risk Score


# model building
fit1 <- coxph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group + riskScore, data = rt1)  # CpM
fit2 <- coxph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, data = rt2)              # BpM
fit3 <- coxph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+riskScore, data = rt3)          # ApM




# fit1
pred1 <- predict(fit1, type = "risk")
c1 <- concordance.index(x = pred1, surv.time = rt1$futime, surv.event = rt1$fustat, method = "noether")
print(c1$c.index)
print(c1$lower)
print(c1$upper)

# fit2
pred2 <- predict(fit2, type = "risk")
c2 <- concordance.index(x = pred2, surv.time = rt2$futime, surv.event = rt2$fustat, method = "noether")
print(c2$c.index)
print(c2$lower)
print(c2$upper)
# fit3
pred3 <- predict(fit3, type = "risk")
c3 <- concordance.index(x = pred3, surv.time = rt3$futime, surv.event = rt3$fustat, method = "noether")
print(c3$c.index)
print(c3$lower)
print(c3$upper)
#Compare the merits of different models using the Likelihood Ratio Test
anova(fit2, fit1, test = "LRT")
anova(fit3, fit1, test = "LRT")
anova(fit2, fit3, test = "LRT")
AIC(fit1, fit2, fit3)


#Draw Nomograma Plot

nom1=regplot(fit1,
             clickable=F,
             title="",
             points=T,#Display the score line
             droplines=F,#Show the connection line
             observation=NULL,
             rank=NULL,
             failtime = c(1,3,5),
             showP = F,
             prfail = F) 

nom2=regplot(fit2,
             clickable=F,
             title="",
             points=TRUE,
             droplines=T,
             observation=NULL,
             rank=NULL,
             failtime = c(1,3,5),
             showP = F,
             prfail = F) 

nom3=regplot(fit3,
             clickable=F,
             title="",
             points=T,
             droplines=T,
             observation=NULL,
             rank=NULL,
             failtime = c(1,3,5),
             showP = F,
             prfail = F) 

#####Merge calibration curves####
####CpMcalibration####
#1 year
pdf(file="CpMcalibration.pdf", width=6, height=6)
fit1 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group + riskScore, x=T, y=T, surv=T, data=rt1, time.inc=1)
cal <- calibrate(fit1, cmethod="KM", method="boot", u=1, m=28,B=1000)##  m = 28 indicates that 20% of the total sample size.
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="Firebrick2", sub=F)
#3 year
fit1 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group + riskScore, x=T, y=T, surv=T, data=rt1, time.inc=3)
cal <- calibrate(fit1, cmethod="KM", method="boot", u=3, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="MediumSeaGreen", sub=F, add=T)
#5 year
fit1 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group + riskScore, x=T, y=T, surv=T, data=rt1, time.inc=5)
cal <- calibrate(fit1, cmethod="KM", method="boot", u=5, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("Firebrick3","MediumSeaGreen","NavyBlue"), lwd=3, bty = 'n')
dev.off()

##BpMcalibration####
pdf(file="BpMcalibration.pdf", width=6, height=6)
#1 year
fit2 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, x=T, y=T, surv=T, data=rt2, time.inc=1)
cal <- calibrate(fit2, cmethod="KM", method="boot", u=1, m=28,B=1000)##   m = 28 indicates that 20% of the total sample size.
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="Firebrick2", sub=F)
#3 year
fit2 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender + Group, x=T, y=T, surv=T, data=rt2, time.inc=3)
cal <- calibrate(fit2, cmethod="KM", method="boot", u=3, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="MediumSeaGreen", sub=F, add=T)
#5 year
fit2 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender + Group, x=T, y=T, surv=T, data=rt2, time.inc=5)
cal <- calibrate(fit2, cmethod="KM", method="boot", u=5, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("Firebrick3","MediumSeaGreen","NavyBlue"), lwd=3, bty = 'n')
dev.off()

####ApMcalibration####
pdf(file="ApMcalibration.pdf", width=6, height=6)
#1 year
fit3 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+riskScore, x=T, y=T, surv=T, data=rt3, time.inc=1)
cal <- calibrate(fit3, cmethod="KM", method="boot", u=1, m=28,B=1000)## m = 28 indicates that 20% of the total sample size.
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="Firebrick2", sub=F)
#3 year
fit3 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+riskScore, x=T, y=T, surv=T, data=rt3, time.inc=3)
cal <- calibrate(fit3, cmethod="KM", method="boot", u=3, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="MediumSeaGreen", sub=F, add=T)
#5 year
fit3 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+riskScore, x=T, y=T, surv=T, data=rt3, time.inc=5)
cal <- calibrate(fit3, cmethod="KM", method="boot", u=5, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("Firebrick3","MediumSeaGreen","NavyBlue"), lwd=3, bty = 'n')
dev.off()



#Read the risk input file
rt=read.table("all.txt", header=T, sep="\t", check.names=F, row.names=1)

dd <- datadist(rt)
options(datadist = "dd")
#risk$Age=as.numeric(risk$Age)

#Read the risk input file
colnames(rt)
rt1=rt[,-c(6,7,8,9,11)]#metastasis、Age、Gender and riskScore
rt2=rt1[,-6]#no riskScore
rt3 <- rt1[,-5]#No metastasis information


# model building
fit1 <- coxph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group + riskScore, data = rt1)  
fit2 <- coxph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, data = rt2)              
fit3 <- coxph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+riskScore, data = rt3)          


#calculate C-index（Harrell's concordance index）

# fit1
pred1 <- predict(fit1, type = "risk")
c1 <- concordance.index(x = pred1, surv.time = rt1$futime, surv.event = rt1$fustat, method = "noether")
print(c1$c.index)
print(c1$lower)
print(c1$upper)

# fit2
pred2 <- predict(fit2, type = "risk")
c2 <- concordance.index(x = pred2, surv.time = rt2$futime, surv.event = rt2$fustat, method = "noether")
print(c2$c.index)
print(c2$lower)
print(c2$upper)
# fit3
pred3 <- predict(fit3, type = "risk")
c3 <- concordance.index(x = pred3, surv.time = rt3$futime, surv.event = rt3$fustat, method = "noether")
print(c3$c.index)
print(c3$lower)
print(c3$upper)
#Compare the advantages and disadvantages of the models using the Likelihood Ratio Test
anova(fit2, fit1, test = "LRT")
anova(fit3, fit1, test = "LRT")
anova(fit2, fit3, test = "LRT")
AIC(fit1, fit2, fit3)

