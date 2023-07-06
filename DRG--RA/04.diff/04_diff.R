

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("reshape2")
#install.packages("ggpubr")


#???ð?
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

expFile="drgGeneExp.txt"      #?????????ļ?
setwd("C:\\biowolf\\geoCRG\\07.diff")      #???ù???Ŀ¼

#??ȡ?????????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
exp=data

#??ȡ??Ʒ?ķ?????Ϣ
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))

#????????????
sigVec=c()
sigGeneVec=c()
for(i in row.names(data)){
	test=wilcox.test(data[i,] ~ Type)
	pvalue=test$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	if(pvalue<0.05){
	sigVec=c(sigVec, paste0(i, Sig))
	sigGeneVec=c(sigGeneVec, i)}
}
#?????????????ı???��
data=data[sigGeneVec,]
outTab=rbind(ID=colnames(data), data)
write.table(outTab, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)
row.names(data)=sigVec

#?Բ??????????п??ӻ?????????ͼ
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=7, height=4.5)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()

#?ѱ???????ת????ggplot2?????ļ?
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=Type)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#????????ͼ
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
	     xlab="",
	     ylab="Gene expression",
	     legend.title="Type",
	     palette = c("blue", "red"),
	     add="point",
	     width=0.8)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#????????ͼ
pdf(file="boxplot.pdf", width=7, height=5)
print(p1)
dev.off()


####