#script deal with the BJ cohort
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(DGEobj.utils)
library(RColorBrewer)

############################################################################
################### read datasets (two batches) ############################
############################################################################
setwd("D:/06_CRC/bulkRNAseq/BeiJingSample/")
GeneCount_Batch1=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/rawData/BeiJing_FirstBatch.CountMatrix.txt",sep="\t",header=T,row.names=1)
GeneCount_Batch2=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/rawData/BeiJing_SecondBatch.CountMatrix.txt",sep="\t",header=T,row.names=1)
GeneList=intersect(rownames(GeneCount_Batch1),rownames(GeneCount_Batch2))
GeneInfo=read.csv("D:/06_CRC/bulkRNAseq/BeiJingSample/GeneInformation.txt",sep="\t",header=T,row.names=1)
GeneList=intersect(GeneList,rownames(GeneInfo))
length(GeneList)
# 61644

GeneCount_Batch1=GeneCount_Batch1[GeneList,]
GeneCount_Batch2=GeneCount_Batch2[GeneList,]
GeneInfo=GeneInfo[GeneList,]

all(rownames(GeneCount_Batch1)==rownames(GeneCount_Batch2))
#TRUE
GeneCount=cbind(GeneCount_Batch1,GeneCount_Batch2)
dim(GeneCount)
#[1] 61644   270

############################################################################
################### read datasets assocaited information ################## 
############################################################################
SampleInformation_raw=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/rawData/SampleInfo_Raw.txt",sep="\t",header=T,row.names=1)
dim(SampleInformation_raw)
#[1] 288   9
### Step 1: Remove samples with no positive or negative information: 4 samples (LN)
SampleInformation_raw=SampleInformation_raw[!is.na(SampleInformation_raw$LNGroup),]
#284   9
### Step 1: Remove samples not sequenced
SampleList=intersect(rownames(SampleInformation_raw),colnames(GeneCount))
length(SampleList)
#266

############################################################################
################### obtain the normalizaed gene expression################## 
############################################################################
GeneCount=GeneCount[,SampleList]
dim(GeneCount)
#61644   266
GeneTPM=convertCounts(
  as.matrix(GeneCount),
  unit="TPM",
  geneLength=GeneInfo$gene_length,
  log = FALSE,
  normalize = "none",
  prior.count = NULL
)

############################################################################
################### sample and gene filter ################## 
############################################################################
### Step 2: Remove samples with distince gene expression
logTPM=log2(GeneTPM+1)
#sampleDist=dist(t(logTPM))
sampleDist=as.dist(1-cor(logTPM))
sampleTree=hclust(sampleDist)
plot(sampleTree,hang=-1)
sampleLable=cutree(sampleTree,k=7)
table(sampleLable)
#  1   2   3   4   5   6   7 
#252  13   1   1   1   1   1
OutLierSample=names(sampleLable[sampleLable>2])
#"LN_317"  "LN2_317" "L3K53"   "L4D6"    "T_317"
KeepSample=names(sampleLable[sampleLable%in%c(1,2)])
length(KeepSample)
#261

GeneTPM_1stFilter=GeneTPM[,KeepSample]
SampleInformation_1stFilter=SampleInformation_raw[KeepSample,]
dim(SampleInformation_1stFilter)
#261   9

### Step 3: Remove LN samples with distince gene expression
LNInfo=SampleInformation_1stFilter[SampleInformation_1stFilter$LNGroup%in%c("Pos","Neg"),]
dim(LNInfo)
#176   9
LNTPM=GeneTPM[,rownames(LNInfo)]
dim(LNTPM)
#61644   176

logLNTPM=log2(LNTPM+1)
#sampleDist=dist(t(logTPM))
LNSampleDist=as.dist(1-cor(logLNTPM))
LNSampleTree=hclust(LNSampleDist)
plot(LNSampleTree,hang=-1)
LNSampleLable=cutree(LNSampleTree,k=2)
table(LNSampleLable)
OutLierLNSample=names(LNSampleLable[LNSampleLable==2])
#"LN_392" "LN_396" "K62L1" 
KeepLNSample=names(LNSampleLable[LNSampleLable==1])
length(KeepLNSample)
#173

### generate the final LN TPM and Count matrix 
LN_TPM_Final=GeneTPM[,KeepLNSample]
dim(LN_TPM_Final)
#61644   173

LN_Count_Final=GeneCount[,KeepLNSample]
dim(LN_Count_Final)
#61644   173

write.table(LN_TPM_Final,file="D:/06_CRC/bulkRNAseq/BeiJingSample/LN_TPM_Final_Ensemble.txt",sep="\t",quote=F)
write.table(LN_Count_Final,file="D:/06_CRC/bulkRNAseq/BeiJingSample/LN_Count_Final_Ensemble.txt",sep="\t",quote=F)

############################################################################
################### tissue sample and gene filter ################## 
############################################################################
TissueInfo=SampleInformation_1stFilter[SampleInformation_1stFilter$LNGroup%in%c("Tissue"),]
dim(TissueInfo)
#85  9
TissueTPM=GeneTPM[,rownames(TissueInfo)]
dim(TissueTPM)
#61644    85
### Step 2: Remove LN samples with distince gene expression
logTissueTPM=log2(TissueTPM+1)
#sampleDist=dist(t(logTPM))
TissueSampleDist=dist(1-cor(logTissueTPM))
TissueSampleTree=hclust(TissueSampleDist)
plot(TissueSampleTree,hang=-1)
KeepTissueSample=colnames(TissueTPM)

KeepSample=c(KeepLNSample,KeepTissueSample)
length(KeepSample)
#258
GeneTPM_Final=GeneTPM[,KeepSample]
GeneCount_Final=GeneCount[,KeepSample]

write.table(GeneTPM_Final,file="D:/06_CRC/bulkRNAseq/BeiJingSample/AllSample_Gene_TPM_Final_Ensemble.txt",sep="\t",quote=F)
write.table(GeneCount_Final,file="D:/06_CRC/bulkRNAseq/BeiJingSample/AllSample_Gene_Count_Final_Ensemble.txt",sep="\t",quote=F)

SampleInformation=SampleInformation_raw[KeepSample,]
write.table(SampleInformation,file="D:/06_CRC/bulkRNAseq/BeiJingSample/SampleInformation_Final.txt",sep="\t",quote=F)

########################################################################################################
######## remove duplicated gene terms with multiple gene_name ##########################
########################################################################################################
GeneTPM_Final_Ensemble=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/AllSample_Gene_TPM_Final_Ensemble.txt",header=T,row.names=1,sep="\t")
GeneInfo=GeneInfo[!duplicated(GeneInfo$gene_name),]
GeneList=intersect(rownames(GeneTPM_Final_Ensemble),rownames(GeneInfo))
length(GeneList)
#58885
GeneInfo=GeneInfo[GeneList,]
GeneTPM_Final_Symbol=GeneTPM_Final_Ensemble[GeneList,]
all(rownames(GeneTPM_Final_Symbol)==rownames(GeneInfo))
rownames(GeneTPM_Final_Symbol)=GeneInfo$gene_name
write.table(GeneTPM_Final_Symbol,file="D:/06_CRC/bulkRNAseq/BeiJingSample/AllSample_Gene_TPM_Final_Symbol.txt",sep="\t",quote=F)

write.table(GeneTPM_Final_Symbol[,KeepLNSample],file="D:/06_CRC/bulkRNAseq/BeiJingSample/LN_TPM_Final_Symbol.txt",sep="\t",quote=F)


########################################################################################################
#########  PCA analyis of all samples (inlucidng lymph nodes and Tissue)  ##########################
########################################################################################################
GeneTPM_Final_Symbol=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/AllSample_Gene_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
SampleInformation=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SampleInformation_Final.txt",sep="\t",header=T,row.names=1)
logTPM=log2(GeneTPM_Final_Symbol+1)
logTPM=logTPM[rowMeans(logTPM)>0,]
logTPM.pca <- prcomp(t(logTPM))
head(logTPM.pca$x)
all(rownames(logTPM.pca$x)==rownames(SampleInformation))

SampleInformation$GroupTmp=paste0(SampleInformation$LNGroup,"_",SampleInformation$Tissue,sep="")
table(SampleInformation$GroupTmp)
#Neg_LN        Pos_LN Tissue_Normal  Tissue_Tumor 
#   139            34            21            64
SampleInformation[SampleInformation$GroupTmp=="Tissue_Normal","Group"]="Normal tissue"
SampleInformation[SampleInformation$GroupTmp=="Tissue_Tumor","Group"]="Tumor tissue"
SampleInformation[SampleInformation$GroupTmp=="Neg_LN","Group"]="Non-metastatic LNs"
SampleInformation[SampleInformation$GroupTmp=="Pos_LN","Group"]="Metastatic LNs"
SampleInformation$Group=factor(SampleInformation$Group,levels=c("Non-metastatic LNs","Metastatic LNs","Normal tissue","Tumor tissue"))
table(SampleInformation$Group)

df_pca=data.frame(PC1=logTPM.pca$x[,1],PC2=logTPM.pca$x[,2],Group=SampleInformation$Group)
#df_pca$Status=ifelse(SampleInformation$Tissue=="LN","Tumor",SampleInformation$Tissue)
table(df_pca$Group)
g=ggplot(df_pca,aes(x=PC1,y=PC2,color=Group,shape=Group))+ 
geom_point(size=3)+
theme_bw()+
scale_shape_manual(values=c(19, 8,5,9))+
scale_color_manual(values=c('SpringGreen4',"Orange",'Khaki', 'Firebrick3'))
pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/AllSampe_PCA_Point.pdf",width=7,height=4.5)
print(g)
dev.off()

eigs <- logTPM.pca$sdev^2
eigs[1] / sum(eigs) #0.3681
eigs[2] / sum(eigs) #0.1471

summary(logTPM.pca)

########################################################################################################
#########  DEGs between PLN and NLN  ##########################
########################################################################################################
LN_Count=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/LN_Count_Final_Ensemble.txt",sep="\t",header=T,row.names=1)

GeneInfo=read.csv("D:/06_CRC/bulkRNAseq/BeiJingSample/GeneInformation.txt",sep="\t",header=T,row.names=1)
GeneInfo=GeneInfo[!duplicated(GeneInfo$gene_name),]
GeneList=intersect(rownames(LN_Count),rownames(GeneInfo))
length(GeneList)
#58885
GeneInfo=GeneInfo[GeneList,]
LN_Count=LN_Count[GeneList,]
all(rownames(LN_Count)==rownames(GeneInfo))
rownames(LN_Count)=GeneInfo$gene_name
dim(LN_Count)
#58884   173
write.table(LN_Count,file="D:/06_CRC/bulkRNAseq/BeiJingSample/LN_Count_Final_Symbol.txt",sep="\t",quote=F)

LN_Count=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/LN_Count_Final_Symbol.txt",header=T,row.names=1)
SampleInformation=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SampleInformation_Final.txt",sep="\t",header=T,row.names=1)
LN_Info=SampleInformation[colnames(LN_Count),]
all(rownames(LN_Info)==colnames(LN_Count))
#TRUE
table(LN_Info$LNGroup)
#Neg Pos 
#139  34
LN_Info$LNGroup=factor(LN_Info$LNGroup,levels=c("Neg","Pos"))
PosvsNeg_dds <- DESeqDataSetFromMatrix(LN_Count, LN_Info, design = ~LNGroup)
PosvsNeg_dds <- DESeq(PosvsNeg_dds)
PosvsNeg_dds <- results(PosvsNeg_dds)
PosvsNeg_dds=PosvsNeg_dds[order(PosvsNeg_dds$padj),]
PosvsNeg_dds=na.omit(PosvsNeg_dds)
PosvsNeg_dds_SigGene=PosvsNeg_dds[PosvsNeg_dds$padj<0.01&abs(PosvsNeg_dds$log2FoldChange)>1,]
table(PosvsNeg_dds_SigGene$Pattern)
write.table(PosvsNeg_dds,file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_PosvsNeg.txt",sep="\t",quote=F)

PosvsNeg_dds=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_PosvsNeg.txt",header=T,row.name=1,sep="\t")
library(ggrepel)
hs_data=data.frame(PosvsNeg_dds)
hs_data$threshold = as.factor(ifelse(hs_data$padj>0.01,"NoSig",ifelse(abs(hs_data$log2FoldChange)>log2(2),ifelse(hs_data$log2FoldChange>log2(2),"SigUp","SigDown"),ifelse(hs_data$log2FoldChange>0,"Up","Down"))))
table(hs_data$threshold)
#Down   NoSig SigDown   SigUp      Up 
#4420   21760    1871    2959    3539
hs_data$threshold=factor(hs_data$threshold,levels=c("SigDown","Down","NoSig","Up","SigUp"))
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = log2FoldChange, y = -log10(padj), size=log2(baseMean+1),colour=threshold, label =ID )) +
  geom_point(alpha=0.4) +
  theme_bw() + scale_size(range = c(0, 4))+
  scale_color_manual(values=c("blue","RoyalBlue", "grey","Salmon","red")) +
  #xlim(c(-4, 4)) + ylim(0,25)+
  geom_vline(xintercept=c(-log2(2),log2(2)),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
    geom_text_repel(
    data = subset(hs_data, hs_data$padj < 0.01 & abs(hs_data$log2FoldChange) >= log2(2)),
    aes(label = ""),
    size = 3,
    max.overlaps=3,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_PosvsNeg_Volcano_NoLabel.pdf",width=9)
print(t)
dev.off()

write.table(hs_data[hs_data$threshold%in%c("SigUp","SigDown"),],file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_PosvsNeg_SigGene.txt",quote=F,sep="\t")
write.table(rownames(hs_data[hs_data$threshold=="SigUp",]),file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_PosvsNeg_SigUp.txt",row.names=F,col.names=F,quote=F)
write.table(rownames(hs_data[hs_data$threshold=="SigDown",]),file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_PosvsNeg_SigDown.txt",row.names=F,col.names=F,quote=F)


LN_TPM=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/LN_TPM_Final_Symbol.txt",header=T,row.names=1)
SampleInformation=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SampleInformation_Final.txt",sep="\t",header=T,row.names=1)
LN_Info=SampleInformation[colnames(LN_TPM),]
all(rownames(LN_Info)==colnames(LN_TPM))
#TRUE
PosvsNeg_dds=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_PosvsNeg.txt",header=T,row.names=1)
PosvsNeg_dds_SigGene=PosvsNeg_dds[PosvsNeg_dds$padj<0.01&abs(PosvsNeg_dds$log2FoldChange)>1,]
PosvsNeg_dds_SigGene$Pattern=ifelse(PosvsNeg_dds_SigGene$log2FoldChange>0,"Up","Down")
table(PosvsNeg_dds_SigGene$Pattern)
pheatmap(LN_TPM[rownames(PosvsNeg_dds_SigGene),],scale="row",show_rownames=FALSE,annotation_row=PosvsNeg_dds_SigGene[,c("log2FoldChange","Pattern")],annotation_col=LN_Info[,c("LNGroup","Batch")])


##############################################################################
#########  SubType identification of Lymphnode  ##########################
##############################################################################
LN_TPM=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/LN_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
SampleInformation=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SampleInformation_Final.txt",sep="\t",header=T,row.names=1)

LN_Info=SampleInformation[colnames(LN_TPM),]
LN_TPM=LN_TPM[rowMeans(LN_TPM)>0,]
LN_TPM_log=log2(LN_TPM+1)
LNSampleDist=as.dist(1-cor(LN_TPM_log))

#SampleInformation_top=read.csv("D:/06_CRC/bulkRNAseq/BeiJingSample/LymphNodeInfomation-Beijing.txt",header=T,row.names=1,sep="\t")
#pheatmap(as.matrix(LNSampleDist),clustering_method="ward.D2",annotation_col=SampleInformation_top[,c("LNMStatus","GroupByGene","Group","Batch")],show_rownames=F,show_colnames=F)

t=pheatmap(as.matrix(LNSampleDist),clustering_method="ward.D2",annotation_col=LN_Info[,c("LNGroup","Batch")],show_rownames=F,show_colnames=F)
SampleAnno=data.frame(Cluster=cutree(t$tree_col,k=7))
all(rownames(SampleAnno)==rownames(LN_Info))
LN_Info$Group=paste0("C",SampleAnno$Cluster,sep="")
pheatmap(as.matrix(LNSampleDist),clustering_method="ward.D2",annotation_col=LN_Info[,c("Group","LNGroup","Batch")],show_rownames=F,show_colnames=F)

LN_Info$GroupByGene=NA
LN_Info[LN_Info$Group=="C7","GroupByGene"]="NLN_C1"
LN_Info[LN_Info$Group=="C1","GroupByGene"]="NLN_C2"
LN_Info[LN_Info$Group=="C4","GroupByGene"]="NLN_C3"
LN_Info[LN_Info$Group%in%c("C2","C6"),"GroupByGene"]="NLN_C4"
LN_Info[LN_Info$Group=="C5","GroupByGene"]="PLN_C1"
LN_Info[LN_Info$Group=="C3","GroupByGene"]="PLN_C2"

GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info$GroupByGene))
ann_colors = list(
                    GroupByGene = GroupByGene_Colors, 
                    Batch = c("Batch1" = "white", Batch2 = "Gray"),
                    LNGroup=c("Pos"="Orange","Neg"="LightBlue")
                )

pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/SubType_Heatmap.pdf",width=7,height=5)
pheatmap(as.matrix(LNSampleDist),clustering_method="ward.D2",
  annotation_colors = ann_colors,annotation_col=LN_Info[,c("GroupByGene","LNGroup","Batch")],
  show_rownames=F,show_colnames=F,color = colorRampPalette(brewer.pal(n = 7, name ="PiYG"))(100))
dev.off()

pheatmap(as.matrix(LNSampleDist),clustering_method="ward.D2",
  annotation_colors = ann_colors,annotation_col=LN_Info[,c("GroupByGene","LNGroup","Batch")],
  show_rownames=F,show_colnames=F,color = colorRampPalette(c("red","firebrick3","lightblue","white","white","white"))(50))

############################# vilidation using various methods #####################################################
SampleDistByEuclidean=dist(t(LN_TPM_log),method = "euclidean")
pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/SubType_Heatmap_euclidean_ward.D2.pdf",width=7,height=5)
pheatmap(as.matrix(SampleDistByEuclidean),clustering_method="ward.D2",
  annotation_colors = ann_colors,annotation_col=LN_Info[,c("GroupByGene","LNGroup","Batch")],
  show_rownames=F,show_colnames=F,color = colorRampPalette(brewer.pal(n = 7, name ="PiYG"))(100))
dev.off()

SampleDistByManhattan=dist(t(LN_TPM_log),method = "manhattan")
pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/SubType_Heatmap_manhattan_ward.D2.pdf",width=7,height=5)
pheatmap(as.matrix(SampleDistByEuclidean),clustering_method="ward.D2",
  annotation_colors = ann_colors,annotation_col=LN_Info[,c("GroupByGene","LNGroup","Batch")],
  show_rownames=F,show_colnames=F,color = colorRampPalette(brewer.pal(n = 7, name ="PiYG"))(100))
dev.off()

SampleDistByPearson=as.dist(1-cor(LN_TPM_log))
pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/SubType_Heatmap_Pearson_Complete.pdf",width=7,height=5)
pheatmap(as.matrix(SampleDistByPearson),clustering_method="complete",
  annotation_colors = ann_colors,annotation_col=LN_Info[,c("GroupByGene","LNGroup","Batch")],
  show_rownames=F,show_colnames=F,color = colorRampPalette(brewer.pal(n = 7, name ="PiYG"))(100))
dev.off()

SampleDistByPearson=as.dist(1-cor(LN_TPM_log))
pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/SubType_Heatmap_Pearson_Average.pdf",width=7,height=5)
pheatmap(as.matrix(LNSampleDist),clustering_method="average",
  annotation_colors = ann_colors,annotation_col=LN_Info[,c("GroupByGene","LNGroup","Batch")],
  show_rownames=F,show_colnames=F,color = colorRampPalette(brewer.pal(n = 7, name ="PiYG"))(100))
dev.off()

LN_Info=LN_Info[order(LN_Info$GroupByGene),]
write.table(LN_Info,file="D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/LN_SubType_Info.txt",sep="\t",quote=F)

########################################################################################################
###################################  PCA analyis of  lymph nodes  ##########################
########################################################################################################
LN_TPM=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/LN_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
LN_Info=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/LN_SubType_Info.txt",header=T,row.names=1,sep="\t")
LN_TPM=LN_TPM[,rownames(LN_Info)]

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

pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/LN_SubType_PCA_Point.pdf",width=6,height=4)
print(g)
dev.off()

eigs <- logTPM.pca$sdev^2
eigs[1] / sum(eigs) #0.3449098
eigs[2] / sum(eigs) #0.1182276
summary(logTPM.pca)

##############################################################################
########################  DEGs between Subtypes of NLN  ######################
##############################################################################
LN_Count=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/LN_Count_Final_Symbol.txt",header=T,row.names=1)
LN_Info=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/LN_SubType_Info.txt",sep="\t",header=T,row.names=1)
LN_Count=LN_Count[,rownames(LN_Info)]
all(rownames(LN_Info)==colnames(LN_Count))
#TRUE
table(LN_Info$GroupByGene)
#NLN_C1 NLN_C2 NLN_C3 NLN_C4 PLN_C1 PLN_C2 
#    50     35     31     32     11     14

GroupByGeneList=as.character(sort(unique(LN_Info$GroupByGene)))
LN_Info$GroupByGene=factor(LN_Info$GroupByGene,levels=sort(unique(LN_Info$GroupByGene)))
PosvsNeg_dds <- DESeqDataSetFromMatrix(LN_Count, LN_Info, design = ~GroupByGene)
PosvsNeg_dds <- DESeq(PosvsNeg_dds)
colData(PosvsNeg_dds)
resultsNames(PosvsNeg_dds)

for(LNSubGroup in GroupByGeneList[-1]){
   LNSubGroupvsC1_dds <- results(PosvsNeg_dds,contrast=c("GroupByGene",LNSubGroup, "NLN_C1"))
   #LNSubGroupvsC1_dds <- results(PosvsNeg_dds,name="GroupByGene_NLN_C2_vs_NLN_C1")
   LNSubGroupvsC1_dds=LNSubGroupvsC1_dds[order(LNSubGroupvsC1_dds$padj),]
   LNSubGroupvsC1_dds=na.omit(LNSubGroupvsC1_dds)
   write.table(LNSubGroupvsC1_dds,file=paste0("D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_",LNSubGroup,"vsNLN_C1.txt",sep=""),sep="\t",quote=F)
}

##############################################################################
########################  Visualization of DEGs from NLN subtypes  ######################
##############################################################################
DEGAll=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_PosvsNeg.txt",header=T,row.names=1)
DEGAll$Gene=rownames(DEGAll)
DEGAll$Group="PLNvsNLN"
GroupByGeneList=c("NLN_C1","NLN_C2","NLN_C3","NLN_C4")
for(LNSubGroup in GroupByGeneList[-1]){
   DEGTmp=read.table(paste0("D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_",LNSubGroup,"vsNLN_C1.txt",sep=""),header=T,row.names=1,sep="\t")
   DEGTmp$Gene=rownames(DEGTmp)
   DEGTmp$Group=paste0(LNSubGroup,"vsC1",sep="")
   DEGAll=rbind(DEGTmp,DEGAll)
}
DEGAll.Sig=DEGAll[DEGAll$padj<0.01&abs(DEGAll$log2FoldChange)>1,]
DEGAll.Sig$Pattern=ifelse(DEGAll.Sig$log2FoldChange>0,"Up","Down")
table(DEGAll.Sig$Group,DEGAll.Sig$Pattern)

t=ggplot(data = DEGAll.Sig, aes(x = Group, y = log2FoldChange, size=-log10(padj),colour=Pattern)) +
  geom_jitter() +
  ylim(-10,10)+
  theme_bw() + scale_size(range = c(0, 3))+
  scale_color_manual(values=c("blue","red")) +
  theme(axis.title.x=element_blank())
pdf("DEG/DEGStatistic_JitterPlot.pdf",width=7,height=5)
print(t)
dev.off()

DEGAll.Sig=DEGAll.Sig[order(DEGAll.Sig$Group,DEGAll.Sig$Pattern,DEGAll.Sig$padj),]
write.table(DEGAll.Sig,file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_sigDEG.txt",sep="\t",quote=F)

########### overlap DEGs between NLN subtypes  ###########
library("ggvenn")
UpGene=DEGAll.Sig[DEGAll.Sig$Pattern=="Up",]
UpGeneList=split(DEGAll.Sig$Gene,DEGAll.Sig$Group)
pdf("DEG/DEGStatistic_Up_venn.pdf",width=6,height=6)
ggvenn(UpGeneList,fill_color = colorRampPalette(brewer.pal(8, "Oranges"))(4))
dev.off()
DownGene=DEGAll.Sig[DEGAll.Sig$Pattern=="Down",]
DownGeneList=split(DEGAll.Sig$Gene,DEGAll.Sig$Group)
pdf("DEG/DEGStatistic_Down_venn.pdf",width=6,height=6)
ggvenn(DownGeneList,fill_color = colorRampPalette(brewer.pal(8, "Blues"))(4))
dev.off()


########### overlap gene expression trends among NLN subtypes  ###########
library(DEGreport)
packageVersion("DEGreport")#1.34.0
library(RColorBrewer)

LN_TPM=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/LN_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
LN_Info=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/LN_SubType_Info.txt",header=T,row.names=1,sep="\t")
LN_TPM=LN_TPM[,rownames(LN_Info)]

NLN_C4vsNLN_C1_DEG=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_NLN_C4vsNLN_C1.txt",header=T,row.names=1,sep="\t")
NLN_C4vsNLN_C1_SigDEG=NLN_C4vsNLN_C1_DEG[NLN_C4vsNLN_C1_DEG$padj<0.01&abs(NLN_C4vsNLN_C1_DEG$log2FoldChange)>1,]
NLN_C4vsNLN_C1_SigDEG$Pattern=ifelse(NLN_C4vsNLN_C1_SigDEG$log2FoldChange>0,"Up","Down")
table(NLN_C4vsNLN_C1_SigDEG$Pattern)
#Down   Up 
#7737 7523
UpGeneList=rownames(NLN_C4vsNLN_C1_SigDEG[NLN_C4vsNLN_C1_SigDEG$Pattern=="Up",])
UpGeneExpr=LN_TPM[UpGeneList,]
#7523  173
all(colnames(UpGeneExpr)==rownames(LN_Info))
#TRUE
LN_Info$GroupByGene=factor(LN_Info$GroupByGene,levels=sort(unique(LN_Info$GroupByGene)))
UpPattern <- degPatterns(UpGeneExpr, LN_Info, time = "GroupByGene",minc = 15)
write.table(UpPattern$normalized,file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEGPattern/NLN_C4vsNLN_C1_UpGenePattern.txt",sep="\t",quote=F)
UpPatternTable=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/DEGPattern/NLN_C4vsNLN_C1_UpGenePattern.txt",header=T)
pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/DEGPattern/NLN_C4vsNLN_C1_UpGenePattern.pdf",height=6,width=10)
degPlotCluster(UpPatternTable,time="GroupByGene",color="GroupByGene",boxes = FALSE,min_genes =  1)+geom_line(aes_string(group="genes"),alpha=0.5)+scale_color_brewer(palette="YlOrRd")+theme_bw()
dev.off()

DownGeneList=rownames(NLN_C4vsNLN_C1_SigDEG[NLN_C4vsNLN_C1_SigDEG$Pattern=="Down",])
DownGeneExpr=LN_TPM[DownGeneList,]
#7523  173
all(colnames(DownGeneExpr)==rownames(LN_Info))
#TRUE
LN_Info$GroupByGene=factor(LN_Info$GroupByGene,levels=sort(unique(LN_Info$GroupByGene)))
DownPattern <- degPatterns(DownGeneExpr, LN_Info, time = "GroupByGene",minc = 15)
write.table(DownPattern$normalized,file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEGPattern/NLN_C4vsNLN_C1_DownGenePattern.txt",sep="\t",quote=F)
DownPatternTable=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/DEGPattern/NLN_C4vsNLN_C1_DownGenePattern.txt",header=T)
pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/DEGPattern/NLN_C4vsNLN_C1_DownGenePattern.pdf",height=6,width=8)
degPlotCluster(DownPatternTable,time="GroupByGene",color="GroupByGene",boxes = FALSE,min_genes =  1)+geom_line(aes_string(group="genes"),alpha=0.5)+scale_color_brewer(palette="YlGn")+theme_bw()
dev.off()


##############################################################################
#####################  Visualization of target genes ##########################
##############################################################################
setwd("D:/06_CRC/bulkRNAseq/BeiJingSample/")
LN_TPM_BeiJing=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/LN_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
LN_Info_BeiJing=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/LN_SubType_Info.txt",header=T,row.names=1,sep="\t")
LN_Info_BeiJing=LN_Info_BeiJing[colnames(LN_TPM_BeiJing),]
all(rownames(LN_Info_BeiJing)==colnames(LN_TPM_BeiJing))
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info_BeiJing$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info_BeiJing$GroupByGene))

targetGenes=c("CD4","MS4A1")
GraphData=cbind(LN_Info_BeiJing[,c("LNGroup","GroupByGene")],t(LN_TPM_BeiJing[targetGenes,]))
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
GraphData=cbind(LN_Info_BeiJing[,c("LNGroup","GroupByGene")],t(LN_TPM_BeiJing[targetGenes,]))
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
GraphData=cbind(LN_Info_BeiJing[,c("LNGroup","GroupByGene")],t(LN_TPM_BeiJing[targetGenes,]))
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
pdf("D:/06_CRC/Graph/Part4/SubType.marker.boxplot.pdf",width=8,height=3)
print(g)
dev.off()

targetGenes=c("ACKR1","KDR","RBP7","PROX1")
GraphData=cbind(LN_Info_BeiJing[,c("LNGroup","GroupByGene")],t(LN_TPM_BeiJing[targetGenes,]))
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
GraphData=cbind(LN_Info_BeiJing[,c("LNGroup","GroupByGene")],t(LN_TPM_BeiJing[targetGenes,]))
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
pdf("D:/06_CRC/Graph/Part4/TF.boxplot.BeiJing.pdf",width=8,height=4)
print(g)
dev.off()

targetGenes=c("CDKN1A","CDKN2A","CDKN2B")
GraphData=cbind(LN_Info_BeiJing[,c("LNGroup","GroupByGene")],t(LN_TPM_BeiJing[targetGenes,]))
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
pdf("D:/06_CRC/Graph/Part4/Senescence.boxplot.BeiJing.pdf",width=3.5,height=5)
print(g)
dev.off()

Marker=read.table("D:/06_CRC/LNscRNAseq/CellTypeIdentification/LN.CellType.Marker.txt",header=T,sep="\t")
SenescenceMarkers=read.table("D:/06_CRC/bulkRNAseq/Senescence/55markersBySegura.txt",sep="\t",header=T)
SenescenceMarkers=SenescenceMarkers[!SenescenceMarkers$GeneNames%in%unique(Marker$gene),]

targetGenes=intersect(SenescenceMarkers[SenescenceMarkers$Group=="Up",]$GeneNames,rownames(LN_TPM_BeiJing))
GraphData=cbind(LN_Info_BeiJing[,c("LNGroup","GroupByGene")],t(LN_TPM_BeiJing[targetGenes,]))
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
