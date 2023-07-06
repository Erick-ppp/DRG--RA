

#install.packages("rms")
#install.packages("rmda")


#???ð?
library(rms)
library(rmda)

inputFile="normalize.txt"             #?????????ļ?
geneFile="importanceGene.XGB.txt"     #?????б??ļ?
setwd("D:\\桌面\\机器学习\\158geoDRG\\09.nomo\\23.Nomo")      #???ù???Ŀ¼

#??ȡ?????ļ?
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(data)=gsub("-", "_", row.names(data))

#??ȡ?????б??ļ?,??ȡ?????????????ı???��
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#??ȡ??Ʒ??????Ϣ
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)
paste(colnames(data), collapse="+")

#???ݴ???
ddist=datadist(rt)
options(datadist="ddist")

#????ģ?ͣ?????????ͼ
lrmModel=lrm(Type~ ENC1IQGAP1+MYH9+CAPZB+FLNB+MYL6ta=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
	fun.at=c(0.0001,0.1,0.3,0.5,0.7,0.9,0.99),
	lp=F, funlabel="Risk of Disease")
#????????ͼ
pdf("Nomo.pdf", width=8, height=6)
p4ot(nomo)
dev.off()

#????У׼????
cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5.5, height=5.5)
plot(cali,
	xlab="Predicted probability",
	ylab="Actual probability", sub=F)
dev.off()

#???ƾ???????
rt$Type=ifelse(rt$Type=="Control", 0, 1)
dc=decision_curve(Type ~ ENC1IQGAP1+MYH9+CAPZB+FLNB+MYL6ta=rt, 
	family = binomial(link ='logit'),
	thresholds= seq(0,1,by = 0.01),
	confidence.intervals = 0.95)
#????DCAͼ??
pdf(file="DCA.pdf", width=5.5, height=5.5)
plot_decision_curve(dc,
	curve.names="Model",
	xlab="Threshold probability",
	cost.benefit.axis=T,
	col="red",
	confidence.intervals=FALSE,
	standardize=FALSE)
dev.off()


####