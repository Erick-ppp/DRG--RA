
#install.packages("caret")
#install.packages("DALEX")
#install.packages("ggplot2")
#install.packages("randomForest")
#install.packages("kernlab")
#install.packages("pROC")
#???ð?
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)



#???ð?
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(pROC)
library(xgboost)
set.seed(123)      #????????
inputFile="diffGeneExp.txt"      #?????ļ?
setwd("D:\\桌面\\机器学习\\12.model")      #???ù???Ŀ¼

#??ȡ?????ļ?
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group

#????ɭ????ģ??
control=trainControl(method="repeatedcv", number = 5, savePredictions=TRUE)
mod_rf = train(Type ~ .,data = data, method='rf', trControl = control)

#????ѧϰģ??
mod_svm=train(Type ~., data = data, method = "svmRadial", prob.model = TRUE, trControl=control)

#XGBģ??
mod_xgb=train(Type ~., data = data, method = "xgbDART", trControl=control,   verbosity = 0)

#GLMģ??
mod_glm=train(Type ~., data = data, method = "glm", family="binomial", trControl=control)

#????Ԥ?⺯??
p_fun=function(object, newdata){
	predict(object, newdata=newdata, type="prob")[,2]
}
yTest=ifelse(data$Type=="con", 0, 1)

#????ɭ????ģ??Ԥ??????
explainer_rf=explain(mod_rf, label = "RF",
                         data = data, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_rf=model_performance(explainer_rf)
#????ѧϰģ??Ԥ??????
explainer_svm=explain(mod_svm, label = "SVM",
                         data = data, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_svm=model_performance(explainer_svm)

#XGBģ??Ԥ??????
explainer_xgb=explain(mod_xgb, label = "XGB",
                      data = data, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_xgb=model_performance(explainer_xgb)
#GLMģ??Ԥ??????
explainer_glm=explain(mod_glm, label = "GLM",
                      data = data, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_glm=model_performance(explainer_glm)

#???Ʋв??ķ????ۼƷֲ?ͼ
pdf(file="residual.pdf", width=6, height=6)
p1 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm)
print(p1)
dev.off()

#???Ʋв???????ͼ
pdf(file="boxplot.pdf", width=6, height=6)
p2 <- plot(mp_rf, mp_svm,  mp_xgb, mp_glm, geom = "boxplot")
print(p2)
dev.off()


#????ROC????
pred1=predict(mod_rf, newx=data, type="prob")
pred2=predict(mod_svm, newx=data, type="prob")
pred3=predict(mod_xgb, newdata=data, type="prob")
pred4=predict(mod_glm, newdata=data, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
roc2=roc(yTest, as.numeric(pred2[,2]))
roc3=roc(yTest, as.numeric(pred3[,2]))
roc4=roc(yTest, as.numeric(pred4[,2]))
ci1=ci.auc(roc1, method="bootstrap")
ci2=ci.auc(roc2, method="bootstrap")
ci3=ci.auc(roc3, method="bootstrap")
ci4=ci.auc(roc4, method="bootstrap")
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="red")
plot(roc2, print.auc=F, legacy.axes=T, main="", col="blue", add=T)
plot(roc3, print.auc=F, legacy.axes=T, main="", col="green", add=T)
plot(roc4, print.auc=F, legacy.axes=T, main="", col="yellow", add=T)
legend('bottomright',
	   c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
	     paste0('SVM: ',sprintf("%.03f",roc2$auc)),
	     paste0('XGB: ',sprintf("%.03f",roc3$auc)),
	     paste0('GLM: ',sprintf("%.03f",roc4$auc))),
	   col=c("red","blue","green","yellow"), lwd=2, bty = 'n')


#?????ַ??????л???????Ҫ?Է???,?õ????ַ?????????Ҫ???��?
importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)
importance_xgb<-variable_importance(
  explainer_xgb,
  loss_function = loss_root_mean_square
)
#???ƻ?????Ҫ??ͼ??
pdf(file="importance.pdf", width=7, height=10)
plot(importance_rf[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_svm[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_xgb[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_glm[c(1,(ncol(data)-8):(ncol(data)+1)),])
dev.off()
#??????Ҫ???��????ߵĻ???
geneNum=5     #???û???????Ŀ
write.table(importance_rf[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.RF.txt", sep="\t", quote=F, row.names=F)
write.table(importance_svm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.SVM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_xgb[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.XGB.txt", sep="\t", quote=F, row.names=F)
write.table(importance_glm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.GLM.txt", sep="\t", quote=F, row.names=F)


dev.off()


