#The relationship between lymph node size and subtype should be clearly demarcated as its own analysis, rather than being inferred from broader trends.
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

BeiJingInfo=read.table("D:/06_CRC/bulkRNAseq/LNInfo/LymphNodeInfomation-Beijing_173.txt",header=T,sep="\t")
BeiJingInfo$CombineGroup=ifelse(BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2"),"NLN_C1/C2",ifelse(BeiJingInfo$GroupByGene%in%c("NLN_C3","NLN_C4"),"NLN_C3/C4","PLN"))

BeiJingInfo_SubType=BeiJingInfo[BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2","NLN_C3","NLN_C4"),TargetInfo]
GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246","#E2BF8C","#147286")
names(GroupByGene_Colors)=sort(unique(BeiJingInfo$GroupByGene))

MSS.MSIStatus=data.frame(table(BeiJingInfo$GroupByGene,BeiJingInfo$MSS.MSI))
colnames(MSS.MSIStatus)=c("SubTypes","MSS.MSI","Number")
MSS.MSIStatus$MSS.MSI=factor(MSS.MSIStatus$MSS.MSI,levels=c("Small","Large"))
g=ggplot(MSS.MSIStatus,aes(x=MSS.MSI, y=Number, fill=SubTypes)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_GroupByGene_MSS.MSI.barplot.pdf",width=3,height=3.5)
print(g)
dev.off()
MSS.MSIStatus.mt=as.data.frame.matrix(table(BeiJingInfo_SubType$GroupByGene,BeiJingInfo_SubType$MSS.MSI))
fisher.test(MSS.MSIStatus.mt) #p-value = 0.07643



##############################################################################
#########################            BeiJing         #########################
##############################################################################
MajorTypeList=c("CD4 T cell","CD8 T cell","NK cell","Exhausted T cell","B cell","Germinal center B cell","Plasma cell","MKI67 cell","Mast cell","Neutrophial","Macrophage","Epithelial cell","Fibroblast cell","Endothelial cell")
BeiJing_LN_Info=read.table("D:/06_CRC/bulkRNAseq/LNInfo/LymphNodeInfomation-Beijing_173.txt",header=T,sep="\t")
rownames(BeiJing_LN_Info)=BeiJing_LN_Info$Sample.ID
BeiJing_Cibersoftx=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/Cibersoftx/CIBERSORTx_Job2_Adjusted.txt",header=T,row.names=1,sep="\t",check.names=F)
#BeiJing_LN_Info=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/LN_SubType_Info.txt",header=T,row.names=1,sep="\t",check.names=F)
LNList=intersect(rownames(BeiJing_Cibersoftx),rownames(BeiJing_LN_Info))
length(LNList)
BeiJing_LN_Info=BeiJing_LN_Info[LNList,]
BeiJing_Cibersoftx=BeiJing_Cibersoftx[LNList,]
all(rownames(BeiJing_Cibersoftx)==rownames(BeiJing_LN_Info))
BeiJing_Cibersoftx=BeiJing_Cibersoftx[,c(1:(ncol(BeiJing_Cibersoftx)-3))]

GraphData.df=cbind(Sample=rownames(BeiJing_LN_Info),MSS.MSI=BeiJing_LN_Info$MSS.MSI,GroupByGene=BeiJing_LN_Info$GroupByGene,BeiJing_Cibersoftx)
GraphData.ldf=reshape2::melt(GraphData.df,id=c(1:3))
colnames(GraphData.ldf)=c("Sample","MSS.MSI","GroupByGene","CellType","CellRatio")

g=ggplot(GraphData.ldf,aes(x=MSS.MSI, y=CellRatio, fill=MSS.MSI)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, alpha=0.7,size=2)+
      facet_wrap(.~CellType,scales="free_y",ncol=5)+
      scale_fill_manual(values=c("Sienna1","Firebrick4")) +
      theme_bw()+
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
      stat_compare_means(comparisons =list(c("pMMR","dMMR")))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/Graph/Part2/BeiJing_Cibersoftx_MMR.boxplot.pdf",width=10,height=5.5)
print(g)
dev.off()




########## DEG between dMMR and pMMR #######################
setwd("D:/06_CRC/bulkRNAseq/BeiJingSample/")
LN_TPM_BeiJing=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/LN_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
LN_Info_BeiJing=read.table("D:/06_CRC/bulkRNAseq/LNInfo/LymphNodeInfomation-Beijing_173.txt",header=T,sep="\t")
rownames(LN_Info_BeiJing)=LN_Info_BeiJing$Sample.ID
LN_Info_BeiJing=LN_Info_BeiJing[colnames(LN_TPM_BeiJing),]
all(rownames(LN_Info_BeiJing)==colnames(LN_TPM_BeiJing))

targetGenes=c("CR1","CR2","FCAMR","FDCSP","MS4A1","VCAN","COL1A1","PDPN")
GraphData=cbind(LN_Info_BeiJing[,c("GroupByGene","MSS.MSI")],t(LN_TPM_BeiJing[targetGenes,]))
GraphData.df=reshape2::melt(GraphData,id=c(1:2))
colnames(GraphData.df)=c("GroupByGene","MSS.MSI","GeneSymbol","GeneExpr")
g=ggplot(GraphData.df,aes(x=MSS.MSI, y=log2(GeneExpr+1), fill=MSS.MSI)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, alpha=0.7,size=2)+
      facet_wrap(.~GeneSymbol,scales="free_y",ncol=5)+
      scale_fill_manual(values=c("Sienna1","Firebrick4")) +
      theme_bw()+
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
      stat_compare_means(comparisons =list(c("pMMR","dMMR")))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/Graph/Part4/SubType.FDCMarker.MSSStatus.boxplot.pdf",width=10,height=4.5)
print(g)
dev.off()





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













