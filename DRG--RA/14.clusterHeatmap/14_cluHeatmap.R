

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

clusterFile="cluster.txt"      #???͵Ľ????ļ?
setwd("C:\\biowolf\\geoCRG\\15.clusterHeatmap")      #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[order(rt$Cluster),]

#׼????ͼ????
data=t(rt[,1:(ncol(rt)-1),drop=F])
Type=rt[,ncol(rt),drop=F]

#??????ͼע?͵???ɫ
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
crgCluCol=bioCol[1:length(levels(factor(Type$Cluster)))]
names(crgCluCol)=levels(factor(Type$Cluster))
ann_colors[["Cluster"]]=crgCluCol

#??????ͼ
pdf("heatmap.pdf", width=7, height=4.5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
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

#??????ת????ggplot2?????ļ?
data=melt(rt, id.vars=c("Cluster"))
colnames(data)=c("Cluster", "Gene", "Expression")

#????????ͼ
p=ggboxplot(data, x="Gene", y="Expression", color = "Cluster",
	     xlab="",
	     ylab="Gene expression",
	     legend.title="Cluster",
	     palette = crgCluCol,
	     width=0.8,
	     add="point")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Cluster),
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#????????ͼ
pdf(file="boxplot.pdf", width=7, height=5)
print(p1)
dev.off()


