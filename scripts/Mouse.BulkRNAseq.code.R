library(ggplot2)
library(RColorBrewer)
library(DEGreport)
MouseData.count=read.csv("D:/06_CRC/bulkRNAseq/Mouse/rawData/gene_count_matrix.csv",row.names=1,header=T)
write.table(colnames(MouseData.count),file="D:/06_CRC/bulkRNAseq/Mouse/rawData/SampleID.txt",sep="\t",row.names=F,quote=F)
dim(MouseData.count)
#57132    71

MouseData.tpm=read.csv("D:/06_CRC/bulkRNAseq/Mouse/rawData/gene_tpm_matrix.csv",row.names=1,header=T)
dim(MouseData.tpm)
#57132    70

GeneNames=data.frame(do.call(rbind, strsplit(rownames(MouseData.count), "\\|")))
MouseData.count$Symbol=GeneNames[,2]
MouseData.count=MouseData.count[!duplicated(MouseData.count$Symbol),]
#56897    71
rownames(MouseData.count)=MouseData.count$Symbol
MouseData.count$Symbol=NULL
#56897    70

GeneNames=data.frame(do.call(rbind, strsplit(rownames(MouseData.tpm), "\\|")))
MouseData.tpm$Symbol=GeneNames[,2]
MouseData.tpm=MouseData.tpm[!duplicated(MouseData.tpm$Symbol),]
dim(MouseData.tpm)
#56897    71
rownames(MouseData.tpm)=MouseData.tpm$Symbol
MouseData.tpm$Symbol=NULL
#56897    70
MouseData.tpm=MouseData.tpm[rowMeans(MouseData.tpm)>0,]

###Remove samples with distince gene expression
logTPM=log2(MouseData.tpm+1)
#sampleDist=dist(t(logTPM))
sampleDist=as.dist(1-cor(logTPM))
sampleTree=hclust(sampleDist)
plot(sampleTree,hang=-1)
# EG6_MNL1

SampleInfo=read.table("D:/06_CRC/bulkRNAseq/Mouse/rawData/SampleInfo.txt",header=T,row.names=1,sep="\t")
table(SampleInfo$Treatment,SampleInfo$Time_c)
SampleInfo=SampleInfo[setdiff(rownames(SampleInfo),c("EG6_MLN1")),]
MouseData.tpm.Filter=MouseData.tpm[rowMeans(MouseData.tpm)>0,rownames(SampleInfo)]
dim(MouseData.tpm.Filter)
#41160    69
all(colnames(MouseData.tpm.Filter)==rownames(SampleInfo))
MouseData.tpm.log=log2(MouseData.tpm.Filter+1)
SampleDist=as.dist(1-cor(MouseData.tpm.log))
pheatmap(as.matrix(SampleDist),annotation_col=SampleInfo)


LNInfo=SampleInfo[SampleInfo$Group=="LN",]
LN_TPM=MouseData.tpm[,rownames(LNInfo)]
dim(LN_TPM)
#51757   162
LN_TPM_log=log2(LN_TPM+1)
LNSampleDist=as.dist(1-cor(LN_TPM_log))
pheatmap(as.matrix(LNSampleDist),clustering_method="ward.D2",annotation_col=LNInfo)
#no difference

##https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
##https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/08b_time_course_analyses.html
##https://bioconductor.org/packages/3.7/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
LN_Count=MouseData.count[,rownames(LNInfo)]
all(rownames(LNInfo)==colnames(LN_Count))
LNInfo$Time_c=factor(LNInfo$Time_c,levels=c("2week","4week","6week","8week","10week","12week"))
LNInfo$Treatment=factor(LNInfo$Treatment,levels=c("Control","AOM"))
dds <- DESeqDataSetFromMatrix(countData = LN_Count, colData = LNInfo, design = ~ Treatment + Time_c + Treatment:Time_c)
dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~Treatment + Time_c)
resultsNames(dds_lrt_time)
res_dds <- results(dds_lrt_time)
res_dds_df=data.frame(res_dds)
res_dds_df=na.omit(res_dds_df)
res_dds_sig=res_dds_df[res_dds_df$pvalue<0.01,]
res_dds_sig$Pattern=ifelse(res_dds_sig$log2FoldChange>0,"Up","Down")
table(res_dds_sig$Pattern)
#Down   Up 
#376  462
write.table(res_dds_sig[order(res_dds_sig$Pattern),],file="D:/06_CRC/bulkRNAseq/Mouse/DEGBetweenAOMAndCtrl.txt",sep="\t",quote=F)
t=rbind(toString(rownames(res_dds_sig[res_dds_sig$Pattern=="Down",])),toString(rownames(res_dds_sig[res_dds_sig$Pattern=="Up",])))
write.table(t,file="D:/06_CRC/bulkRNAseq/Mouse/DEGBetweenAOMAndCtrl.GeneList.txt",sep="\t",quote=F,col.names=F)
write.table(rownames(res_dds_sig[res_dds_sig$Pattern=="Down",]),file="D:/06_CRC/bulkRNAseq/Mouse/DownAOM.GeneList.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(rownames(res_dds_sig[res_dds_sig$Pattern=="Up",]),file="D:/06_CRC/bulkRNAseq/Mouse/UpAOM.GeneList.txt",sep="\t",quote=F,col.names=F,row.names=F)


res_dds_12 <- results(dds_lrt_time,name="TreatmentAOM.Time_c12week", test="Wald")
res_dds_df_12=data.frame(res_dds_12)
res_dds_df_12=na.omit(res_dds_df_12)
res_dds_sig_12=res_dds_df_12[abs(res_dds_df_12$log2FoldChange)>log2(1.5),]
dim(res_dds_sig_12)
GeneList=intersect(rownames(res_dds_sig),rownames(res_dds_sig_12))
length(GeneList)

vsd <- vst(dds, blind=FALSE)
GeneExpr=assay(vsd)
LNSampleDist=as.dist(1-cor(GeneExpr))
pheatmap(as.matrix(LNSampleDist),clustering_method="ward.D2",annotation_col=LNInfo[,c("Time","Group","Treatment")])
#no difference
pheatmap(GeneExpr[rownames(res_dds_sig),],clustering_method="ward.D2",annotation_col=LNInfo[,c("Time","Group","Treatment")],scale="row")






UpGeneExpr=GeneExpr[rownames(res_dds_sig[res_dds_sig$Pattern=="Up",]),]
dim(UpGeneExpr)
LNInfo$Time_c=factor(LNInfo$Time_c,levels=c("2week","4week","6week","8week","10week","12week"))
UpPattern <- degPatterns(UpGeneExpr, metadata = LNInfo, time="Time_c", col="Treatment",minc =  15)
write.table(UpPattern$normalized,file="D:/06_CRC/bulkRNAseq/Mouse/UpDEGPattern.txt",sep="\t",quote=F)
UpPatternTable=read.table("D:/06_CRC/bulkRNAseq/Mouse/UpDEGPattern.txt",header=T)
UpPatternTable$Time_c=factor(UpPatternTable$Time_c,levels=c("2week","4week","6week","8week","10week","12week"))
degPlotCluster(UpPatternTable,time="Time_c",color="Treatment",min_genes =  1)+
scale_color_manual(values=c("orange","lightgrey"))+theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
pdf("D:/06_CRC/bulkRNAseq/Mouse/UpDEGPattern.pdf",height=6,width=8)
degPlotCluster(UpPatternTable,time="Time_c",color="Treatment",facet =F,lines =F)+
scale_color_manual(values=c("red","lightgrey"))+
theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
dev.off()
UpGene=unique(UpPatternTable[,c("genes","cluster")])
write.table(UpGene[order(UpGene$cluster),],file="D:/06_CRC/bulkRNAseq/Mouse/UpDEGPattern_ClusterInfo.txt",sep="\t",quote=F,row.names=F)



DownGeneExpr=GeneExpr[rownames(res_dds_sig[res_dds_sig$Pattern=="Down",]),]
LNInfo$Time_c=factor(LNInfo$Time_c,levels=c("2week","4week","6week","8week","10week","12week"))
DownPattern <- degPatterns(DownGeneExpr, metadata = LNInfo, time="Time_c", col="Treatment",minc =  1)
write.table(DownPattern$normalized,file="D:/06_CRC/bulkRNAseq/Mouse/DownDEGPattern.txt",sep="\t",quote=F)
DownPatternTable=read.table("D:/06_CRC/bulkRNAseq/Mouse/DownDEGPattern.txt",header=T)
DownPatternTable$Time_c=factor(DownPatternTable$Time_c,levels=c("2week","4week","6week","8week","10week","12week"))
degPlotCluster(DownPatternTable,time="Time_c",color="Treatment",min_genes =  10)+
scale_color_manual(values=c("RoyalBlue","lightgrey"))+theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))

DownPatternTable_filter=DownPatternTable[DownPatternTable$cluster%in%c(1,10,16,18,23,5,7,8),]
degPlotCluster(DownPatternTable_filter,time="Time_c",color="Treatment",min_genes =  10)+
scale_color_manual(values=c("RoyalBlue","lightgrey"))+theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))

pdf("D:/06_CRC/bulkRNAseq/Mouse/DownDEGPattern.pdf",height=6,width=8)
degPlotCluster(DownPatternTable_filter,time="Time_c",color="Treatment",facet =F,lines =F)+
scale_color_manual(values=c("RoyalBlue","lightgrey"))+theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))
dev.off()

length(DownGene)
write.table(DownGene[order(DownGene$cluster),],file="D:/06_CRC/bulkRNAseq/Mouse/DownDEGPattern_ClusterInfo.txt",sep="\t",quote=F,row.names=F)

write.table(unique(DownPatternTable_filter$genes),file="D:/06_CRC/bulkRNAseq/Mouse/DownAOM.GeneList.ByPattern.txt",sep="\t",quote=F,row.names=F,col.names=F)
length(unique(DownPatternTable_filter$genes))

targetGenes=rownames(res_dds_sig)[1:16]
targetGenes=c("Adamts1","Agtr1a","Apoe","Bmp2","Cav1","Cdkn1b","Comt","Edn1","Egr1","Esr1","S1pr2","Hbegf","Jun","Klf4","Nos3","Npr3","Ptgs2","Thbs1","Timp3","Tnfaip3","Vegfa","Ndrg2")
targetGenes=c("Adm","Agtr1a","Apoe","Nr2f2","Rhob","Bmp2","Zfp36l1","Cav1","Cxcr4","Comt","Edn1","Efna1","Egr1","Ephb2")
targetGenes=c("Cd4","Cd8a","Ms4a1","Cd79","Col1a1","Acta2","Fn1","Pdpn","Col5a1","Vcan","Fdcsp","Cr1","Cr2","Cdkn1a","Cdkn2a","Cdkn1b")
targetGeneExpr=GeneExpr[intersect(targetGenes,rownames(GeneExpr)),]
all(rownames(LNInfo)==colnames(targetGeneExpr))
targetGeneExprInfo=cbind(LNInfo,t(targetGeneExpr))
targetGeneExprInfo.df=reshape2::melt(targetGeneExprInfo,id=c(1:ncol(LNInfo)))
colnames(targetGeneExprInfo.df)=c(colnames(LNInfo),"Symbol","Expr")
targetGeneExprInfo.df$Time_c=factor(targetGeneExprInfo.df$Time_c,levels=c("2week","4week","6week","8week","10week","12week"))
ggplot(targetGeneExprInfo.df,
  aes(x = Time_c, y = Expr, color = Treatment, group = Treatment)) + 
  facet_wrap(Symbol~.,scales="free_y")+
  scale_color_manual(values=c("lightgrey","Orange"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))+
  geom_point() + geom_smooth(se = FALSE, method = "loess")










TissueInfo=SampleInfo[SampleInfo$Group=="Tissue",]
Tissue_TPM=MouseData.tpm[,rownames(TissueInfo)]
dim(Tissue_TPM)
Tissue_Count=MouseData.count[,rownames(TissueInfo)]
all(rownames(TissueInfo)==colnames(Tissue_Count))
TissueInfo$Time_c=factor(TissueInfo$Time_c,levels=c("2week","4week","6week","8week","10week","12week"))
TissueInfo$Treatment=factor(TissueInfo$Treatment,levels=c("Control","AOM"))
dds <- DESeqDataSetFromMatrix(countData = Tissue_Count, colData = TissueInfo, design = ~ Treatment + Time_c + Treatment:Time_c)
dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~Treatment + Time_c)
resultsNames(dds_lrt_time)
res_dds <- results(dds_lrt_time)
res_dds_df=data.frame(res_dds)
res_dds_df=na.omit(res_dds_df)
res_dds_sig=res_dds_df[res_dds_df$pvalue<0.01,]
res_dds_sig$Pattern=ifelse(res_dds_sig$log2FoldChange>0,"Up","Down")
table(res_dds_sig$Pattern)
#Down   Up 
#376  462
write.table(res_dds_sig[order(res_dds_sig$Pattern),],file="D:/06_CRC/bulkRNAseq/Mouse/Tisse_DEGBetweenAOMAndCtrl.txt",sep="\t",quote=F)
write.table(rownames(res_dds_sig[res_dds_sig$Pattern=="Down",]),file="D:/06_CRC/bulkRNAseq/Mouse/Tisse_DownAOM.GeneList.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(rownames(res_dds_sig[res_dds_sig$Pattern=="Up",]),file="D:/06_CRC/bulkRNAseq/Mouse/Tisse_UpAOM.GeneList.txt",sep="\t",quote=F,col.names=F,row.names=F)

Tissue_vsd <- vst(dds, blind=FALSE)
Tissue_GeneExpr=assay(Tissue_vsd)

targetGenes=rownames(res_dds_sig)[1:16]
targetGeneExpr=Tissue_GeneExpr[intersect(targetGenes,rownames(Tissue_GeneExpr)),]
all(rownames(TissueInfo)==colnames(targetGeneExpr))
targetGeneExprInfo=cbind(TissueInfo,t(targetGeneExpr))
targetGeneExprInfo.df=reshape2::melt(targetGeneExprInfo,id=c(1:ncol(TissueInfo)))
colnames(targetGeneExprInfo.df)=c(colnames(TissueInfo),"Symbol","Expr")
targetGeneExprInfo.df$Time_c=factor(targetGeneExprInfo.df$Time_c,levels=c("2week","4week","6week","8week","10week","12week"))
ggplot(targetGeneExprInfo.df,
  aes(x = Time_c, y = Expr, color = Treatment, group = Treatment)) + 
  facet_wrap(Symbol~.,scales="free_y")+
  scale_color_manual(values=c("lightgrey","Orange"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))+
  geom_point() + geom_smooth(se = FALSE, method = "loess")




