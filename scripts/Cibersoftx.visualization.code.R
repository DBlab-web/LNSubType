library(ggplot2)
library(ggpubr)
library(RColorBrewer)

##############################################################################
#########################            BeiJing         #########################
##############################################################################
MajorTypeList=c("CD4 T cell","CD8 T cell","NK cell","Exhausted T cell","B cell","Germinal center B cell","Plasma cell","MKI67 cell","Mast cell","Neutrophial","Macrophage","Epithelial cell","Fibroblast cell","Endothelial cell")


BeiJing_Cibersoftx=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/Cibersoftx/CIBERSORTx_Job2_Adjusted.txt",header=T,row.names=1,sep="\t",check.names=F)
BeiJing_LN_Info=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/LN_SubType_Info.txt",header=T,row.names=1,sep="\t",check.names=F)
LNList=intersect(rownames(BeiJing_Cibersoftx),rownames(BeiJing_LN_Info))
BeiJing_LN_Info=BeiJing_LN_Info[LNList,]
BeiJing_Cibersoftx=BeiJing_Cibersoftx[LNList,]
all(rownames(BeiJing_Cibersoftx)==rownames(BeiJing_LN_Info))
BeiJing_Cibersoftx=BeiJing_Cibersoftx[,c(1:(ncol(BeiJing_Cibersoftx)-3))]

GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(BeiJing_LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(BeiJing_LN_Info$GroupByGene))

GraphData.df=cbind(Sampple=rownames(BeiJing_LN_Info),GroupByGene=BeiJing_LN_Info$GroupByGene,BeiJing_Cibersoftx)
GraphData.ldf=reshape2::melt(GraphData.df,id=c(1:2))
colnames(GraphData.ldf)=c("Sample","GroupByGene","CellType","CellRatio")
GraphData.ldf$CellType=factor(GraphData.ldf$CellType,levels=MajorTypeList)
CellType_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(MajorTypeList))
GraphData.ldf$GroupByGene=factor(GraphData.ldf$GroupByGene,levels=c(sort(unique(GraphData.ldf$GroupByGene))))
g=ggplot(GraphData.ldf,aes(x=GroupByGene, y=CellRatio, fill=GroupByGene)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, size=0.1)+
      facet_wrap(.~CellType,scales="free_y",ncol=8)+
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
      stat_compare_means(comparisons =list(c("NLN_C2","NLN_C1"),c("NLN_C3","NLN_C1"),c("NLN_C4","NLN_C1")))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/Graph/Part2/BeiJing_Cibersoftx_Result.boxplot.pdf",width=20,height=5.5)
print(g)
dev.off()

g=ggplot(GraphData.ldf,aes(x=Sample, y=CellRatio, fill=CellType)) +
      geom_bar(position = 'fill',stat = "identity") +
      facet_grid(.~GroupByGene,scales="free_x",space ="free_x")+
      scale_fill_manual(values=CellType_colors) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_blank())
pdf("D:/06_CRC/Graph/Part2/BeiJing_Cibersoftx_Result.barplot.pdf",width=12,height=2)
print(g)
dev.off()


##############################################################################
#########################           ShanXi           #########################
##############################################################################

ShanXi_Cibersoftx=read.table("D:/06_CRC/bulkRNAseq/ShanXiSample/Cibersoftx/CIBERSORTx_Job3_Adjusted.txt",header=T,row.names=1,sep="\t",check.names=F)
ShanXi_LN_Info=read.table("D:/06_CRC/bulkRNAseq/ShanXiSample/LN_GroupInfo.txt",header=T,row.names=1,sep="\t")
LNList=intersect(rownames(ShanXi_Cibersoftx),rownames(ShanXi_LN_Info))
ShanXi_LN_Info=ShanXi_LN_Info[LNList,]
ShanXi_Cibersoftx=ShanXi_Cibersoftx[LNList,]
all(rownames(ShanXi_Cibersoftx)==rownames(ShanXi_LN_Info))
ShanXi_Cibersoftx=ShanXi_Cibersoftx[,c(1:(ncol(ShanXi_Cibersoftx)-3))]

GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(ShanXi_LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(ShanXi_LN_Info$GroupByGene))

GraphData.df=cbind(Sampple=rownames(ShanXi_LN_Info),GroupByGene=ShanXi_LN_Info$GroupByGene,ShanXi_Cibersoftx)
GraphData.ldf=reshape2::melt(GraphData.df,id=c(1:2))
colnames(GraphData.ldf)=c("Sample","GroupByGene","CellType","CellRatio")
GraphData.ldf$CellType=factor(GraphData.ldf$CellType,levels=MajorTypeList)
CellType_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(GraphData.ldf$CellType)))
GraphData.ldf$GroupByGene=factor(GraphData.ldf$GroupByGene,levels=c(sort(unique(GraphData.ldf$GroupByGene))))
GraphData.ldf=GraphData.ldf[order(GraphData.ldf$GroupByGene),]

g=ggplot(GraphData.ldf,aes(x=GroupByGene, y=CellRatio, fill=GroupByGene)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, size=0.1)+
      facet_wrap(.~CellType,scales="free_y",ncol=8)+
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
      stat_compare_means(comparisons =list(c("NLN_C2","NLN_C1"),c("NLN_C3","NLN_C1"),c("NLN_C4","NLN_C1")))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/Graph/Part2/ShanXi_Cibersoftx_Result.boxplot.pdf",width=20,height=5.5)
print(g)
dev.off()


g=ggplot(GraphData.ldf,aes(x=Sample, y=CellRatio, fill=CellType)) +
      geom_bar(position = 'fill',stat = "identity") +
      facet_grid(.~GroupByGene,scales="free_x",space ="free_x")+
      scale_fill_manual(values=CellType_colors) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_blank())
pdf("D:/06_CRC/Graph/Part2/ShanXi_Cibersoftx_Result.barplot.pdf",width=12,height=2)
print(g)
dev.off()








##############################################################################
#########################           HarBin           #########################
##############################################################################

HarBin_Cibersoftx=read.table("D:/06_CRC/bulkRNAseq/HarBinSample/Cibersoftx/CIBERSORTx_Job4_Adjusted.txt",header=T,row.names=1,sep="\t",check.names=F)
HarBin_LN_Info=read.table("D:/06_CRC/bulkRNAseq/HarBinSample/LN_GroupInfo.txt",header=T,row.names=1,sep="\t")
LNList=intersect(rownames(HarBin_Cibersoftx),rownames(HarBin_LN_Info))
HarBin_LN_Info=HarBin_LN_Info[LNList,]
HarBin_Cibersoftx=HarBin_Cibersoftx[LNList,]
HarBin_Cibersoftx=HarBin_Cibersoftx[,c(1:(ncol(HarBin_Cibersoftx)-3))]
all(rownames(HarBin_Cibersoftx)==rownames(HarBin_LN_Info))

GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(HarBin_LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(HarBin_LN_Info$GroupByGene))

GraphData.df=cbind(Sampple=rownames(HarBin_LN_Info),GroupByGene=HarBin_LN_Info$GroupByGene,HarBin_Cibersoftx)
GraphData.ldf=reshape2::melt(GraphData.df,id=c(1:2))
colnames(GraphData.ldf)=c("Sample","GroupByGene","CellType","CellRatio")
GraphData.ldf$CellType=factor(GraphData.ldf$CellType,levels=MajorTypeList)
CellType_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(GraphData.ldf$CellType)))
GraphData.ldf$GroupByGene=factor(GraphData.ldf$GroupByGene,levels=c(sort(unique(GraphData.ldf$GroupByGene))))
GraphData.ldf=GraphData.ldf[order(GraphData.ldf$GroupByGene),]
g=ggplot(GraphData.ldf,aes(x=GroupByGene, y=CellRatio, fill=GroupByGene)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.25, size=0.1)+
      facet_wrap(.~CellType,scales="free_y",ncol=8)+
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
      stat_compare_means(comparisons =list(c("NLN_C2","NLN_C1"),c("NLN_C3","NLN_C1"),c("NLN_C4","NLN_C1")))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/Graph/Part2/HarBin_Cibersoftx_Result.boxplot.pdf",width=20,height=5.5)
print(g)
dev.off()


g=ggplot(GraphData.ldf,aes(x=Sample, y=CellRatio, fill=CellType)) +
      geom_bar(position = 'fill',stat = "identity") +
      facet_grid(.~GroupByGene,scales="free_x",space ="free_x")+
      scale_fill_manual(values=CellType_colors) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_blank())
pdf("D:/06_CRC/Graph/Part2/HarBin_Cibersoftx_Result.barplot.pdf",width=12,height=2)
print(g)
dev.off()






