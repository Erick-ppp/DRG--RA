

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


library(ConsensusClusterPlus)      #???ð?
expFile="diffGeneExp.txt"          #?????????ļ?
workDir="C:\\biowolf\\geoCRG\\14.cluster"     #????Ŀ¼
setwd(workDir)      #???ù???Ŀ¼

#??ȡ?????ļ?
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

#ȥ??????????Ʒ, ֻ????ʵ??????Ʒ
group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
data=data[,group=="Treat"]

#????
maxK=9     #??????????kֵ
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=123456,
              plot="png")

#һ???Դ???
calcICL(results, title="consensusScore", plot="png")

#???????ͽ???
clusterNum=2        #?ּ??࣬????ǰ????ͼ???ж?
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("Cluster")
cluster$Cluster=paste0("C", cluster$Cluster)
outTab=cbind(t(data), cluster)
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="cluster.txt", sep="\t", quote=F, col.names=F)


