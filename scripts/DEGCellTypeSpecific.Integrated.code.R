#this script was used to estimate the cell type specific gene expression of DEGs between various compairsion, correspondent to Figure 3A
library(presto) #download zip file from github, unzip, install.packages(path,rep=NULL,type="resource")
library(dplyr)
library(msigdbr)
library(tibble)
library(fgsea)

##### LN.Combined.Final.rds was generated from ##############
setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/")
LNCombined.Final=readRDS("LN.Combined.Final.rds")
table(LNCombined.Final$cellType)

MajorType.marker <- wilcoxauc(LNCombined.Final, 'cellType')
table(MajorType.marker$group)

setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample")
DEGAll=read.table("DEG/AllGene_PosvsNeg.txt",header=T,row.names=1)
DEGAll$Gene=rownames(DEGAll)
DEGAll$Group="PLNvsNLN"
GroupByGeneList=c("NLN_C1","NLN_C2","NLN_C3","NLN_C4")
for(LNSubGroup in GroupByGeneList[-1]){
   DEGTmp=read.table(paste0("DEG/AllGene_",LNSubGroup,"vsNLN_C1.txt",sep=""),header=T,row.names=1,sep="\t")
   DEGTmp$Gene=rownames(DEGTmp)
   DEGTmp$Group=paste0(LNSubGroup,"vsC1",sep="")
   DEGAll=rbind(DEGTmp,DEGAll)
}
DEGAll.Sig=DEGAll[DEGAll$padj<0.01&abs(DEGAll$log2FoldChange)>1,]
DEGAll.Sig$Pattern=ifelse(DEGAll.Sig$log2FoldChange>0,"Up","Down")
table(DEGAll.Sig$Group,DEGAll.Sig$Pattern)

DEGAll.Sig=DEGAll.Sig[DEGAll.Sig$Gene%in%rownames(LNCombined.Final),]
DEGAll.Sig$GroupInfo=paste0(DEGAll.Sig$Pattern,"_",DEGAll.Sig$Group,sep="")
fgsea_sets=split(DEGAll.Sig$Gene,DEGAll.Sig$GroupInfo)

for(cluster in unique(LNCombined.Final$cellType)){
print (cluster)
clusterCell<- MajorType.marker %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
ranks=na.omit(ranks)
fgseaRes <- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0, nPermSimple = 10000)
fgseaRes$cellType=cluster
fwrite(fgseaRes, file=paste0("GSEA/DEGSpecificExpr/",cluster,".txt",sep=""), sep="\t", sep2=c("", " ", ""),quote=F)
}


cluster="Fibroblast cell"
clusterCell<- MajorType.marker %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
ranks=na.omit(ranks)
fgseaRes <- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0, nPermSimple = 10000)

pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample/GSEA/Up_NLN_C2vsC1_Fib.pdf",width=3,height=2)
plotEnrichment(fgsea_sets[["Up_NLN_C2vsC1"]],ranks)
dev.off()
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample/GSEA/Up_NLN_C3vsC1_Fib.pdf",width=3,height=2)
plotEnrichment(fgsea_sets[["Up_NLN_C3vsC1"]],ranks)
dev.off()
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample/GSEA/Up_PLNvsNLN_Fib.pdf",width=3,height=2)
plotEnrichment(fgsea_sets[["Up_PLNvsNLN"]],ranks)
dev.off()




cluster="Endothelial cell"
clusterCell<- MajorType.marker %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
ranks=na.omit(ranks)
fgseaRes <- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0, nPermSimple = 10000)
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample/GSEA/Up_NLN_C3vsC1_End.pdf",width=3,height=2)
plotEnrichment(fgsea_sets[["Up_NLN_C3vsC1"]],ranks)
dev.off()
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample/GSEA/Up_PLNvsNLN_End.pdf",width=3,height=2)
plotEnrichment(fgsea_sets[["Up_PLNvsNLN"]],ranks)
dev.off()




MajorTypeList=c("CD4 T cell","CD8 T cell","NK cell","Exhausted T cell","B cell","Germinal center B cell","Plasma cell","MKI67 cell","Mast cell","Neutrophial","Macrophage","Epithelial cell","Fibroblast cell","Endothelial cell")

setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample/GSEA/DEGSpecificExpr")
for(type in MajorTypeList){
	if(type=="CD4 T cell"){
		AllResult=read.table(paste0(type,".txt",sep=""),header=T,sep="\t",check.names=F)
	}else{
		tmp=read.table(paste0(type,".txt",sep=""),header=T,sep="\t",check.names=F)
		AllResult=rbind(AllResult,tmp)
	}
}

AllResult=AllResult[,c("pathway","pval","padj","NES","size","cellType")]
AllResult$cellType=factor(AllResult$cellType,levels=rev(MajorTypeList))
t=na.omit(AllResult$padj)
minp=min(t[!t==0])
AllResult$padj=ifelse(AllResult$padj==0,minp*0.01,AllResult$padj)
t=ggplot(data = AllResult, aes(x = pathway, y = cellType, size=-log10(padj),colour=NES)) +
  geom_point() +
  theme_bw() + 
  scale_size(range = c(1, 3))+
  scale_colour_gradient2(low="white",mid="white",midpoint=0,high="red")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample/GSEA/GSEA_DEGSpecificExpr.pdf",width=5,height=4)
print(t)
dev.off()



