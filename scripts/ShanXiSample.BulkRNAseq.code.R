
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(DESeq2)
library(DGEobj.utils)

GeneCount_Batch1=read.csv("D:/06_CRC/bulkRNAseq/ShanXiSample/rawData/FirstBatch_113Sample/gene_count_matrix.csv",sep=",",header=T,row.names=1)
GeneCount_Batch2=read.csv("D:/06_CRC/bulkRNAseq/ShanXiSample/rawData/SecondBatch_234Sample/gene_count_matrix.csv",sep=",",header=T,row.names=1)
all(rownames(GeneCount_Batch1)==rownames(GeneCount_Batch2))
#FALSE
GeneList=intersect(rownames(GeneCount_Batch1),rownames(GeneCount_Batch2))
length(GeneList)
#60665
GeneCount_Batch1=GeneCount_Batch1[GeneList,]
GeneCount_Batch2=GeneCount_Batch2[GeneList,]
all(rownames(GeneCount_Batch1)==rownames(GeneCount_Batch2))
#TRUE
GeneCount=cbind(GeneCount_Batch1,GeneCount_Batch2)
dim(GeneCount)
#60656   347

GeneNames=data.frame(do.call(rbind, strsplit(rownames(GeneCount), "\\|")))
GeneNames$Ensemble=data.frame(do.call(rbind, strsplit(GeneNames[,1], "\\.")))[,1]

GeneCount$Ensemble=GeneNames$Ensemble
GeneCount$Symbol=GeneNames[,2]
GeneCount=GeneCount[!duplicated(GeneCount$Symbol),]

GeneInfo=read.csv("D:/06_CRC/bulkRNAseq/BeiJingSample/GeneInformation.txt",sep="\t",header=T,row.names=1)
dim(GeneInfo)
#61852     9
GeneInfo=GeneInfo[!duplicated(GeneInfo$gene_name),]
dim(GeneInfo)
# 59086     9
GeneList=intersect(GeneCount$Ensemble,rownames(GeneInfo))
length(GeneList)
#58898

rownames(GeneCount)=GeneCount$Ensemble
GeneCount=GeneCount[GeneList,]
#58898   349
GeneInfo=GeneInfo[GeneList,]

rownames(GeneCount)=GeneCount$Symbol
GeneCount$Symbol=NULL
GeneCount$Ensemble=NULL
#58898   347

all(rownames(GeneCount)==rownames(GeneInfo$gene_name))
##TRUE

GeneTPM=convertCounts(
  as.matrix(GeneCount),
  unit="TPM",
  geneLength=GeneInfo$gene_length,
  log = FALSE,
  normalize = "none",
  prior.count = NULL
)




SampleInfo=read.table("D:/06_CRC/bulkRNAseq/ShanXiSample/LymphNodeInfomationFromShanXiBatch.txt",sep="\t",header=T,row.names=1)
dim(SampleInfo)
# 227   8

### output gene expression matrix
LNList=intersect(rownames(SampleInfo),colnames(GeneTPM))
length(LNList)
#227
LNTPM=GeneTPM[,LNList]
LNCount=GeneCount[,LNList]

###Remove samples with distince gene expression
logTPM=log2(LNTPM+1)
#sampleDist=dist(t(logTPM))
sampleDist=as.dist(1-cor(logTPM))
sampleTree=hclust(sampleDist)
plot(sampleTree,hang=-1)
# did not find any outlier samples

write.table(LNTPM,file="D:/06_CRC/bulkRNAseq/ShanXiSample/LNSample_Gene_TPM_Syumbol.txt",sep="\t",quote=F)
write.table(LNCount,file="D:/06_CRC/bulkRNAseq/ShanXiSample/LNSample_Gene_Count_Syumbol.txt",sep="\t",quote=F)


##############################################################################
#########  SubType validation of Lymphnode  ##########################
##############################################################################
LN_TPM=read.table("D:/06_CRC/bulkRNAseq/ShanXiSample/LNSample_Gene_TPM_Syumbol.txt",header=T,row.names=1,sep="\t")
LN_Info=read.table("D:/06_CRC/bulkRNAseq/ShanXiSample/LymphNodeInfomationFromShanXiBatch.txt",sep="\t",header=T,row.names=1)
LN_Info$LNStatus=ifelse(LN_Info$Group=="NLN","Neg","Pos")

LN_TPM=LN_TPM[rowMeans(LN_TPM)>0,]
dim(LN_TPM)
#51757   162
LN_TPM_log=log2(LN_TPM+1)
LNSampleDist=as.dist(1-cor(LN_TPM_log))

#pheatmap(as.matrix(LNSampleDist),clustering_method="ward.D2",annotation_col=LN_Info[,c("LNStatus","GroupByGene")],show_rownames=F,show_colnames=F)

t=pheatmap(as.matrix(LNSampleDist),clustering_method="ward.D2",annotation_col=LN_Info[,c("LNStatus","GroupByGene")],show_rownames=F,show_colnames=F)
SampleAnno=data.frame(Cluster=cutree(t$tree_col,k=6))
all(rownames(SampleAnno)==rownames(LN_Info))
LN_Info$Group=paste0("C",SampleAnno$Cluster,sep="")
pheatmap(as.matrix(LNSampleDist),clustering_method="ward.D2",annotation_col=LN_Info[,c("Group","LNStatus","GroupByGene")],show_rownames=F,show_colnames=F)

LN_Info$GroupByGene=NA
LN_Info[LN_Info$Group%in%c("C1"),"GroupByGene"]="NLN_C1"
LN_Info[LN_Info$Group=="C2","GroupByGene"]="NLN_C2"
LN_Info[LN_Info$Group=="C3","GroupByGene"]="NLN_C3"
LN_Info[LN_Info$Group=="C4","GroupByGene"]="NLN_C4"
LN_Info[LN_Info$Group%in%c("C5","C6"),"GroupByGene"]="PLN_C1"
LN_Info$LNGroup=LN_Info$LNStatus

GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info$GroupByGene))
ann_colors = list(
                    GroupByGene = GroupByGene_Colors, 
                    LNGroup=c("Pos"="Orange","Neg"="LightBlue")
                )
pdf("D:/06_CRC/bulkRNAseq/ShanXiSample/SubType_Heatmap.pdf",width=7,height=5)
pheatmap(as.matrix(LNSampleDist),clustering_method="ward.D2",
  annotation_colors = ann_colors,annotation_col=LN_Info[,c("GroupByGene","LNGroup")],
  show_rownames=F,show_colnames=F,color = colorRampPalette(brewer.pal(n = 7, name ="PiYG"))(100))
dev.off()

write.table(LN_Info,file="D:/06_CRC/bulkRNAseq/ShanXiSample/LN_GroupInfo.txt",sep="\t",quote=F)


##############################################################################
#########  PCA analyis of both LN and Tissue  ##########################
##############################################################################

LN_TPM=read.table("D:/06_CRC/bulkRNAseq/ShanXiSample/LNSample_Gene_TPM_Syumbol.txt",header=T,row.names=1,sep="\t")
LN_Info=read.table("D:/06_CRC/bulkRNAseq/ShanXiSample/LN_GroupInfo.txt",header=T,row.names=1,sep="\t")
all(rownames(LN_Info)==colnames(LN_TPM))
#TRUE
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info$GroupByGene))

logTPM=log2(LN_TPM+1)
logTPM=logTPM[rowMeans(logTPM)>0,]
logTPM.pca <- prcomp(t(logTPM))
head(logTPM.pca$x)
all(rownames(logTPM.pca$x)==rownames(LN_Info))
df_pca=data.frame(PC1=logTPM.pca$x[,1],PC2=logTPM.pca$x[,2],GroupByGene=LN_Info$GroupByGene,LNStatus=LN_Info$LNGroup)
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info$GroupByGene))
g=ggplot(df_pca,aes(x=PC1,y=PC2,color=GroupByGene,shape=LNStatus))+ 
geom_point(size=3)+
theme_bw()+
scale_shape_manual(values=c(19, 8))+
scale_color_manual(values=GroupByGene_Colors)

pdf("D:/06_CRC/bulkRNAseq/ShanXiSample/ShanXi_LN_SubType_PCA_Point.pdf",width=6,height=4)
print(g)
dev.off()

eigs <- logTPM.pca$sdev^2
eigs[1] / sum(eigs) #0.2917599
eigs[2] / sum(eigs) #0.1128584


summary(logTPM.pca)


##############################################################################
###################################  Visualization  ##########################
##############################################################################

setwd("D:/06_CRC/bulkRNAseq/ShanXiSample/")
LN_TPM_ShanXi=read.table("D:/06_CRC/bulkRNAseq/ShanXiSample/LNSample_Gene_TPM_Syumbol.txt",header=T,row.names=1,sep="\t")
LN_Info_ShanXi=read.table("D:/06_CRC/bulkRNAseq/ShanXiSample/LN_GroupInfo.txt",header=T,row.names=1,sep="\t")
LN_Info_ShanXi=LN_Info_ShanXi[colnames(LN_TPM_ShanXi),]
all(rownames(LN_Info_ShanXi)==colnames(LN_TPM_ShanXi))
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info_ShanXi$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info_ShanXi$GroupByGene))

targetGenes=c("CD4","MS4A1")
GraphData=cbind(LN_Info_ShanXi[,c("LNGroup","GroupByGene")],t(LN_TPM_ShanXi[targetGenes,]))
GraphData.df=reshape2::melt(GraphData,id=c(1:2))
colnames(GraphData.df)=c("LNGroup","GroupByGene","GeneSymbol","GeneExpr")
g=ggplot(GraphData.df,aes(x=GroupByGene, y=log2(GeneExpr+1), fill=GroupByGene)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, alpha=0.6,size=0.5)+
      facet_wrap(.~GeneSymbol,scales="free_y",ncol=1)+
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/Graph/Part2/T.B.marker.boxplot.pdf",width=3.5,height=3)
print(g)
dev.off()

targetGenes=c("COL1A1","ACTA2","FN1","PECAM1","LYVE1")
GraphData=cbind(LN_Info_ShanXi[,c("LNGroup","GroupByGene")],t(LN_TPM_ShanXi[targetGenes,]))
GraphData.df=reshape2::melt(GraphData,id=c(1:2))
colnames(GraphData.df)=c("LNGroup","GroupByGene","GeneSymbol","GeneExpr")
g=ggplot(GraphData.df,aes(x=GroupByGene, y=log2(GeneExpr+1), fill=GroupByGene)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, alpha=0.7,size=0.1)+
      facet_wrap(.~GeneSymbol,scales="free_y",ncol=3)+
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/Graph/Part3/FibAndEnd.marker.boxplot.pdf",width=7,height=3)
print(g)
dev.off()


targetGenes=c("PDPN","COL5A1","VCAN","KDR","FDCSP","CR1","CR2")
GraphData=cbind(LN_Info_ShanXi[,c("LNGroup","GroupByGene")],t(LN_TPM_ShanXi[targetGenes,]))
GraphData.df=reshape2::melt(GraphData,id=c(1:2))
colnames(GraphData.df)=c("LNGroup","GroupByGene","GeneSymbol","GeneExpr")
g=ggplot(GraphData.df,aes(x=GroupByGene, y=log2(GeneExpr+1), fill=GroupByGene)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, alpha=0.7,size=0.1)+
      facet_wrap(.~GeneSymbol,scales="free_y",ncol=4)+
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/Graph/Part4/SubType.marker.boxplot.ShanXi.pdf",width=8,height=3)
print(g)
dev.off()

targetGenes=c("ACKR1","KDR","RBP7","PROX1")
GraphData=cbind(LN_Info_ShanXi[,c("LNGroup","GroupByGene")],t(LN_TPM_ShanXi[targetGenes,]))
GraphData.df=reshape2::melt(GraphData,id=c(1:2))
colnames(GraphData.df)=c("LNGroup","GroupByGene","GeneSymbol","GeneExpr")
ggplot(GraphData.df,aes(x=GroupByGene, y=log2(GeneExpr+1), fill=GroupByGene)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, alpha=0.7,size=0.1)+
      facet_wrap(.~GeneSymbol,scales="free_y",ncol=4)+
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))

targetGenes=c("EGR1","CEBPD","KLF10","ETS2","CEBPB","FOXF2","NFIC","PRNP","SMARCA1","NFIX","SPIB","LTF","RELB","KLF4","ERG","MECOM","SOX18")
GraphData=cbind(LN_Info_ShanXi[,c("LNGroup","GroupByGene")],t(LN_TPM_ShanXi[targetGenes,]))
GraphData.df=reshape2::melt(GraphData,id=c(1:2))
colnames(GraphData.df)=c("LNGroup","GroupByGene","GeneSymbol","GeneExpr")
g=ggplot(GraphData.df,aes(x=GroupByGene, y=log2(GeneExpr+1), fill=GroupByGene)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, alpha=0.7,size=0.1)+
      facet_wrap(.~GeneSymbol,scales="free_y",ncol=5)+
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/Graph/Part4/TF.boxplot.ShanXi.pdf",width=8,height=4)
print(g)
dev.off()


targetGenes=c("CDKN1A","CDKN2A","CDKN2B")
GraphData=cbind(LN_Info_ShanXi[,c("LNGroup","GroupByGene")],t(LN_TPM_ShanXi[targetGenes,]))
GraphData.df=reshape2::melt(GraphData,id=c(1:2))
colnames(GraphData.df)=c("LNGroup","GroupByGene","GeneSymbol","GeneExpr")
g=ggplot(GraphData.df,aes(x=GroupByGene, y=log2(GeneExpr+1), fill=GroupByGene)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, alpha=0.7,size=0.1)+
      facet_wrap(.~GeneSymbol,scales="free_y",ncol=1)+
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/Graph/Part4/Senescence.boxplot.ShanXi.pdf",width=3.5,height=5)
print(g)
dev.off()

