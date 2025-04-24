#The relationship between lymph node size and subtype should be clearly demarcated as its own analysis, rather than being inferred from broader trends.
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

BeiJingInfo=read.table("D:/06_CRC/bulkRNAseq/LNInfo/LymphNodeInfomation-Beijing_173.txt",header=T,sep="\t")
BeiJingInfo$CombineGroup=ifelse(BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2"),"NLN_C1/C2",ifelse(BeiJingInfo$GroupByGene%in%c("NLN_C3","NLN_C4"),"NLN_C3/C4","PLN"))

BeiJingInfo_SubType=BeiJingInfo[BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2","NLN_C3","NLN_C4"),TargetInfo]
GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246","#E2BF8C","#147286")
names(GroupByGene_Colors)=sort(unique(BeiJingInfo$GroupByGene))

LNSizeStatus=data.frame(table(BeiJingInfo$GroupByGene,BeiJingInfo$LNSize))
colnames(LNSizeStatus)=c("SubTypes","LNSize","Number")
LNSizeStatus$LNSize=factor(LNSizeStatus$LNSize,levels=c("Small","Large"))
g=ggplot(LNSizeStatus,aes(x=LNSize, y=Number, fill=SubTypes)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_GroupByGene_LNSize.barplot.pdf",width=3,height=3.5)
print(g)
dev.off()
LNSizeStatus.mt=as.data.frame.matrix(table(BeiJingInfo_SubType$GroupByGene,BeiJingInfo_SubType$LNSize))
fisher.test(LNSizeStatus.mt) #p-value = 0.07643


BeiJingInfo=read.table("D:/06_CRC/bulkRNAseq/LNInfo/LymphNodeInfomation-Beijing_173.txt",header=T,sep="\t")
LN_TPM=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/LN_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
LNList=intersect(colnames(LN_TPM),BeiJingInfo$Sample.ID)
LN_TPM=LN_TPM[,LNList]
LNInfo=BeiJingInfo
rownames(LNInfo)=LNInfo$Sample.ID
LNInfo=LNInfo[LNList,]
all(rownames(LNInfo)==colnames(LN_TPM))

##############################################################################
####################################     All gene   ##########################
##############################################################################
logTPM=log2(LN_TPM+1)
logTPM=logTPM[rowMeans(logTPM)>0,]
logTPM.pca <- prcomp(t(logTPM))
head(logTPM.pca$x)
all(rownames(logTPM.pca$x)==rownames(LNInfo))
df_pca=data.frame(PC1=logTPM.pca$x[,1],PC2=logTPM.pca$x[,2],GroupByGene=LNInfo$GroupByGene,LNSize=LNInfo$LNSize)
#GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LNInfo$GroupByGene)))
GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246","#E2BF8C","#147286")
names(GroupByGene_Colors)=sort(unique(LNInfo$GroupByGene))
g=ggplot(df_pca,aes(x=PC1,y=PC2,color=LNSize,shape=LNSize))+ 
geom_point(size=3)+
theme_bw()+
scale_shape_manual(values=c(19, 8))+
scale_color_manual(values=c("firebrick","LightBlue"))
pdf("D:/06_CRC/Graph/Part6/PCA_LNSize.pdf",width=6,height=4)
      print(g)
dev.off()



##############################################################################
#########################            BeiJing         #########################
##############################################################################
MajorTypeList=c("CD4 T cell","CD8 T cell","NK cell","Exhausted T cell","B cell","Germinal center B cell","Plasma cell","MKI67 cell","Mast cell","Neutrophial","Macrophage","Epithelial cell","Fibroblast cell","Endothelial cell")

BeiJing_Cibersoftx=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/Cibersoftx/CIBERSORTx_Job2_Adjusted.txt",header=T,row.names=1,sep="\t",check.names=F)
#BeiJing_LN_Info=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/LN_SubType_Info.txt",header=T,row.names=1,sep="\t",check.names=F)
BeiJing_LN_Info=LNInfo
LNList=intersect(rownames(BeiJing_Cibersoftx),rownames(BeiJing_LN_Info))
BeiJing_LN_Info=BeiJing_LN_Info[LNList,]
BeiJing_Cibersoftx=BeiJing_Cibersoftx[LNList,]
all(rownames(BeiJing_Cibersoftx)==rownames(BeiJing_LN_Info))
BeiJing_Cibersoftx=BeiJing_Cibersoftx[,c(1:(ncol(BeiJing_Cibersoftx)-3))]

GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246","#E2BF8C","#147286")
names(GroupByGene_Colors)=sort(unique(LNInfo$GroupByGene))
GraphData.df=cbind(Sample=rownames(BeiJing_LN_Info),LNSize=BeiJing_LN_Info$LNSize,GroupByGene=BeiJing_LN_Info$GroupByGene,BeiJing_Cibersoftx)
GraphData.ldf=reshape2::melt(GraphData.df,id=c(1:3))
colnames(GraphData.ldf)=c("Sample","LNSize","GroupByGene","CellType","CellRatio")
GraphData.ldf$CellType=factor(GraphData.ldf$CellType,levels=MajorTypeList)
CellType_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(MajorTypeList))
GraphData.ldf$GroupByGene=factor(GraphData.ldf$GroupByGene,levels=c(sort(unique(GraphData.ldf$GroupByGene))))
GraphData.ldf=na.omit(GraphData.ldf)
g=ggplot(GraphData.ldf[GraphData.ldf$LNSize=="Large",],aes(x=GroupByGene, y=CellRatio, fill=GroupByGene)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, size=0.1)+
      facet_wrap(.~CellType,scales="free_y",ncol=8)+
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
      stat_compare_means(comparisons =list(c("NLN_C2","NLN_C1"),c("NLN_C3","NLN_C1"),c("NLN_C4","NLN_C1")))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/Graph/Part2/BeiJing_Cibersoftx_LNLarge.boxplot.pdf",width=20,height=5.5)
print(g)
dev.off()

g=ggplot(GraphData.ldf[GraphData.ldf$LNSize=="Small",],aes(x=GroupByGene, y=CellRatio, fill=GroupByGene)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, size=0.1)+
      facet_wrap(.~CellType,scales="free_y",ncol=8)+
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
      stat_compare_means(comparisons =list(c("NLN_C2","NLN_C1"),c("NLN_C3","NLN_C1"),c("NLN_C4","NLN_C1")))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/Graph/Part2/BeiJing_Cibersoftx_LNSmall.boxplot.pdf",width=20,height=5.5)
print(g)
dev.off()
