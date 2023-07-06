

#install.packages("caret")
#install.packages("DALEX")
#install.packages("ggplot2")
#install.packages("randomForest")
#install.packages("kernlab")
#install.packages("pROC")
#install.packages("xgboost")


#???ð?
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)

set.seed(123)      #????????
inputFile="test.normalize.txt"         #?????????ļ?
geneFile="importanceGene.XGB.txt"      #?????б??ļ?
setwd("C:\\biowolf\\geoCRG\\25.testROC")      #???ù???Ŀ¼

#??ȡ?????????ļ?
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(data)=gsub("-", "_", row.names(data))

#??ȡ?????б??ļ?, ??ȡ?????????????ı???��
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#??ȡ??Ʒ??????Ϣ
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group

#?????ݽ??з???
inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

#ѡ??ģ??
control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
if(geneFile=="importanceGene.RF.txt"){
	#????ɭ????ģ??
	model=train(Type ~ ., data = train, method='rf', trControl = control)
}else if(geneFile=="importanceGene.SVM.txt"){
	#SVM????ѧϰģ??
	model=train(Type ~., data = train, method = "svmRadial", prob.model=TRUE, trControl=control)
}else if(geneFile=="importanceGene.XGB.txt"){
	#XGBģ??
	model=train(Type ~., data = train, method = "xgbDART", trControl=control)
}else if(geneFile=="importanceGene.GLM.txt"){
	#GLMģ??
	model=train(Type ~., data = train, method = "glm", family="binomial", trControl=control)
}

#????ROC????
yTest=ifelse(test$Type=="Control", 0, 1)
pred1=predict(model, newdata=test, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=T, legacy.axes=T, main="", col="red")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()


####