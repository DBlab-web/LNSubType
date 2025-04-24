###script used to dealt with the scRNA-seq data #####

library(Seurat) #4.4.0
library(harmony) #1.2.1
library(plyr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(scCustomize)

setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq")

PLNSample=readRDS("CRCLD.integrated.cellType.rds")
#21633 features across 84623 samples within 1 assay
NuoHeZhiYinMergeSample=readRDS("merge.rds")
#32245 features across 771675 samples within 1 assay
LNSample=subset(NuoHeZhiYinMergeSample,Type=="LN")
#32245 features across 422136 samples within 1 assay

PLNSample$platform="10XGenomics"
PLNSample$Groups="NA"
PLNSample$Type=ifelse(PLNSample$Sample%in%c("CRC1_T","CRC2_T","CRC3_T"),"Tissue","LN")
merge$Treatment=ifelse(merge$platform=="NuoHe",merge$Group,merge$Groups)

PLNSample@meta.data=PLNSample@meta.data[,c("Sample","nCount_RNA","nFeature_RNA","Treatment","Type","platform","cellType")]

LNSample$Sample=ifelse(is.na(LNSample$original),LNSample$orig.ident,LNSample$original)
LNSample$cellType=LNSample$Celltype
LNSample@meta.data=LNSample@meta.data[,c("Sample","nCount_RNA","nFeature_RNA","Treatment","Type","platform","cellType")]

LNCombined=merge(LNSample,PLNSample)

rm(list=ls())
gc()
LNCombined=readRDS("LN.Combined.raw.rds")
#32270 features across 506759 samples within 1 assay
options(future.globals.maxSize = 1000 * 1024^3)
LNCombined <- NormalizeData(LNCombined, normalization.method = "LogNormalize", scale.factor = 10000)
LNCombined <- FindVariableFeatures(LNCombined, selection.method = "vst", nfeatures = 2000)
LNCombined <- ScaleData(LNCombined)
LNCombined <- RunPCA(LNCombined, features = VariableFeatures(object = LNCombined))
LNCombined <- RunHarmony(LNCombined, group.by.vars=c("platform","Sample"))
LNCombined<- RunUMAP(LNCombined, reduction = "harmony", dims = 1:50)
LNCombined<- FindNeighbors(LNCombined, reduction = "harmony", dims = 1:50)
LNCombined<- FindClusters(LNCombined, resolution = 1)


pdf("CellTypeIdentification/RawMajorCellType.pdf",width=9,height=6)
DimPlot(LNCombined,reduction ="umap",label=TRUE,raster=TRUE,group.by="cellType")&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("CellTypeIdentification/RawCellClusters.pdf",width=9,height=6)
DimPlot(LNCombined,reduction ="umap",label=TRUE,raster=TRUE)&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("CellTypeIdentification/RawPlatform.pdf",width=9,height=6)
DimPlot(LNCombined,reduction ="umap",label=TRUE,group.by="platform",raster=TRUE)&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("CellTypeIdentification/RawSamples.pdf",width=9,height=6)
DimPlot(LNCombined,reduction ="umap",label=TRUE,group.by="Sample",raster=TRUE)&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
saveRDS(LNCombined,file="LN.Combined.rds")


MarkerGene=c("CD3D","CD4","CD8A","CD8B","NKG7","LINC00299","KRT86","MS4A1","LMO2","MEF2B","IGHG3","MZB1","MKI67","MS4A1","KIT","G0S2","LYZ","APOE","EPCAM","CLDN7","COL1A1","ACTA2","CLDN5","ENG")
pdf("CellTypeIdentification/RawCellMarker.pdf",width=20,height=20)
FeaturePlot(LNCombined,features = MarkerGene,raster=TRUE,ncol=6)&NoLegend()&NoAxes()&theme(plot.title=element_text(face="italic"),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()


SampleCluster=data.frame(table(LNCombined$seurat_clusters,LNCombined$Sample))
colnames(SampleCluster)=c("SubCluster","Sample","Number")
ClusterTotal=data.frame(table(LNCombined$seurat_clusters))
colnames(ClusterTotal)=c("SubCluster","TotalNumber")
SampleCluster.df=merge(SampleCluster,ClusterTotal,by="SubCluster")
SampleCluster.df$Ratio=SampleCluster.df$Number/SampleCluster.df$TotalNumber
SampleCluster.df=SampleCluster.df[SampleCluster.df$Ratio>0,]
SampleCluster.df=SampleCluster.df[order(SampleCluster.df$SubCluster,SampleCluster.df$Ratio,decreasing=T),]
SampleCluster.df %>% group_by(SubCluster) %>% top_n(n = 1, wt = Ratio) -> top1
top1Sum <- by(top1$Ratio, top1$SubCluster, sum) #we start by top3, so need sum
top1Sum.df<-data.frame(RatioSum=do.call(rbind,as.list(top1Sum)))
top1Sum.df$Cluster=rownames(top1Sum.df)
SampleBiasSubCluster=top1Sum.df[top1Sum.df$RatioSum>0.5,"Cluster"]
SampleBiasSubCluster
# "15" "30" "32" "33" "35" "36" "37" "38" "39" "40" "41" "42" "43" "44" "45"
# "46" "47" "48" "49" "50" "51" "52" "53" "54" "56" "57" "58"
SampleCluster.df[SampleCluster.df$SubCluster%in%SampleBiasSubCluster,]


### gene expression for each sub cell cluster ################ 
LNCombinedDownsampling=LNCombined[, sample(colnames(LNCombined), size = 10000, replace=F)]
BloodSubCluster.markers <- FindAllMarkers(LNCombinedDownsampling, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(BloodSubCluster.markers,file="CellTypeIdentification/Subcluster.Marker.txt",sep="\t",quote=F)

BloodSubCluster.markers=read.table("CellTypeIdentification/Subcluster.Marker.txt",sep="\t",header=T)
for(i in unique(LNCombined$seurat_clusters)){
t=FeaturePlot(LNCombined, features = BloodSubCluster.markers[BloodSubCluster.markers$cluster==i,"gene"][1:16],raster=TRUE)&NoLegend()&NoAxes()&theme(plot.title=element_text(face="italic"),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
tiff(paste0("CellTypeIdentification/SubSclusterMarker/C",i,"_Marker_Harmony.Subcluster.tiff",sep=""),width=1000,height=1000)
print(t)
dev.off()
}

LNCombined <- PercentageFeatureSet(LNCombined, "^MT-", col.name = "percent_mito")
summary(LNCombined$percent_mito)

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
pdf("CellTypeIdentification/Harmony.nFeaturesRNA.pdf",height=8,width=15)
VlnPlot(LNCombined,group.by = "seurat_clusters", features =feats, pt.size = 0, ncol = 1) +NoLegend()
dev.off()


Idents(LNCombined)=LNCombined$seurat_clusters
new.cluster.ids <- c("B cell","CD4 T cell","B cell","NK cell","Lysing cell","CD4 T cell","CD4 T cell","CD8 T cell","CD4 T cell","CD4 T cell","CD4 T cell","Lysing cell","Plasma cell","Plasma cell","MKI67 cell","Epithelial cell","Macrophage","Germinal center B cell","NK cell","CD4 T cell","CD4 T cell","CD4 T cell","CD8 T cell","Neutrophial","Exhausted T cell","B cell","B cell","Fibroblast cell","Macrophage","Mast cell","Epithelial cell","Endothelial cell","Epithelial cell","CD4 T cell","Doublets","Doublets","Rare","Rare","Rare","Rare","Rare","Rare","Rare","Rare","Rare","Rare","Rare","Rare","Rare","Rare","Rare","Rare","Rare","Rare","Rare","CD4 T cell","Rare","Rare","Rare")
names(new.cluster.ids) <- levels(LNCombined)
LNCombined <- RenameIdents(LNCombined, new.cluster.ids)
LNCombined$cellType=Idents(LNCombined)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(LNCombined$cellType)))
pdf("CellTypeIdentification/LNCombined.Raw.CellType.pdf",width=7,height=5)
DimPlot(LNCombined,cols =mycolors,label=T,raster=TRUE)&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()

saveRDS(LNCombined,file="/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/LNCombined.cellType.rds")

LNCombined.Final=subset(LNCombined,cellType%in%c(setdiff(unique(LNCombined$cellType),c("Lysing cell","Rare","Doublets"))))
LNCombined.Final<-RunUMAP(LNCombined.Final, reduction = "harmony", dims = 1:50)
LNCombined.Final<- FindNeighbors(LNCombined.Final, reduction = "harmony", dims = 1:50)
LNCombined.Final<- FindClusters(LNCombined.Final, resolution = 1)

pdf("CellTypeIdentification/LNCombined.Final.CellType.pdf",width=7,height=5)
DimPlot(LNCombined.Final,label=T,raster=TRUE)&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("CellTypeIdentification/LNCombined.Final.Cluster.pdf",width=7,height=5)
DimPlot(LNCombined.Final,label=T,group.by="seurat_clusters",raster=TRUE)&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
pdf("CellTypeIdentification/Harmony.Final.nFeaturesRNA.pdf",height=8,width=15)
VlnPlot(LNCombined.Final,group.by = "seurat_clusters", features =feats, pt.size = 0, ncol = 1) +NoLegend()
dev.off()
MarkerGene=c("CD3D","CD4","CD8A","CD8B","NKG7","LINC00299","KRT86","MS4A1","LMO2","MEF2B","IGHG3","MZB1","MKI67","MS4A2","KIT","G0S2","LYZ","APOE","EPCAM","CLDN7","COL1A1","ACTA2","CLDN5","ENG")
pdf("CellTypeIdentification/LNCombined.Final.CellMarker.pdf",width=20,height=15)
FeaturePlot(LNCombined.Final,features = MarkerGene,raster=TRUE,ncol=6)&NoLegend()&NoAxes()&theme(plot.title=element_text(face="italic"),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()


new.cluster.ids<-c("B cell","CD4 T cell","B cell","NK cell","CD4 T cell","CD4 T cell","CD4 T cell","CD8 T cell","CD4 T cell","CD4 T cell","Plasma cell","Plasma cell","Epithelial cell","Germinal center B cell","Macrophage","CD4 T cell","NK cell","B cell","MKI67 cell","CD4 T cell","CD8 T cell","CD8 T cell","Neutrophial","Germinal center B cell","Exhausted T cell","B cell","Fibroblast cell","Macrophage","Mast cell","Epithelial cell","Endothelial cell","NK cell","Epithelial cell","CD4 T cell","B cell","CD4 T cell","CD4 T cell","CD4 T cell","CD4 T cell","B cell","CD4 T cell","CD4 T cell")
names(new.cluster.ids) <- levels(LNCombined.Final)
LNCombined.Final <- RenameIdents(LNCombined.Final, new.cluster.ids)
LNCombined.Final$cellType=Idents(LNCombined.Final)
LNCombined.Final$cellType=factor(LNCombined.Final$cellType,levels=c("CD4 T cell","CD8 T cell","NK cell","Exhausted T cell","B cell","Germinal center B cell","Plasma cell","MKI67 cell","Mast cell","Neutrophial","Macrophage","Epithelial cell","Fibroblast cell","Endothelial cell"))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(LNCombined.Final$cellType)))
pdf("CellTypeIdentification/LNCombined.Final.Modified.CellType.pdf",width=6.5,height=4)
DimPlot(LNCombined.Final,cols =mycolors,raster=TRUE,group.by="cellType")&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()


MarkerGene=c("CD3D","CD4","CD8A","CD8B","NKG7","LINC00299","KRT86","MS4A1","LMO2","MEF2B","IGHG3","MZB1","MKI67","MS4A2","KIT","G0S2","LYZ","APOE","EPCAM","CLDN7","COL1A1","ACTA2","CLDN5","ENG")
LNCombined.Final.DownSampling=LNCombined.Final[, sample(colnames(LNCombined.Final), size = 10000, replace=F)]
LNCombined.Final.DownSampling$cellType=factor(LNCombined.Final.DownSampling$cellType,levels=c("CD4 T cell","CD8 T cell","NK cell","Exhausted T cell","B cell","Germinal center B cell","Plasma cell","MKI67 cell","Mast cell","Neutrophial","Macrophage","Epithelial cell","Fibroblast cell","Endothelial cell"))
CellTypeColors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(LNCombined.Final$cellType)))
pdf("CellTypeIdentification/LNCombined.Final.CellType.StackedVlnPlot.pdf",height=15,width=10)
Stacked_VlnPlot(seurat_object = LNCombined.Final.DownSampling,x_lab_rotate=90,plot_spacing =0,group.by="cellType", features = MarkerGene, colors_use = CellTypeColors)
dev.off()

saveRDS(LNCombined.Final,file="LN.Combined.Final.rds")



##############################################################################
############## visualization for each condtions ##############
##############################################################################
setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/")
LNCombined.Final=readRDS("LN.Combined.Final.rds")
#32270 features across 463054 samples within 1 assay
MajorTypeList=c("CD4 T cell","CD8 T cell","NK cell","Exhausted T cell","B cell","Germinal center B cell","Plasma cell","MKI67 cell","Mast cell","Neutrophial","Macrophage","Epithelial cell","Fibroblast cell","Endothelial cell")

LNCombined.Final$LNStatus=ifelse(LNCombined.Final$Sample%in%c("CRC1_MLN1","CRC1_MLN2"),"Pos",ifelse(LNCombined.Final$Sample%in%c("CRC1_T","CRC2_T","CRC3_T"),"Tissue","Neg"))

MajorTypeColors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(LNCombined.Final$cellType)))
names(MajorTypeColors)=MajorTypeList
pdf("CellTypeIdentification/LNCombined.Final.Type.pdf",width=8,height=2.5)
DimPlot(LNCombined.Final,cols =MajorTypeColors,group.by="cellType",raster=TRUE,split.by="LNStatus",ncol=3,pt.size=2)&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()

pdf("CellTypeIdentification/LNCombined.Final.TreatmentNoTissue.pdf",width=8,height=2.5)
DimPlot(subset(LNCombined.Final,Type=="LN"),group.by="cellType",cols =MajorTypeColors,raster=TRUE,split.by="Treatment")&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()


MajorTypeStatus=data.frame(table(LNCombined.Final$LNStatus,LNCombined.Final$cellType))
colnames(MajorTypeStatus)=c("Status","CellType","Number")
MajorTypeStatus$CellType=factor(MajorTypeStatus$CellType,levels=MajorTypeList)
g=ggplot(MajorTypeStatus,aes(x=Status, y=Number, fill=CellType)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=MajorTypeColors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("CellTypeIdentification/LNCombined.Final.Type.CellRatio.pdf",width=4,height=4)
print(g)
dev.off()

MajorTypeStatus=data.frame(table(LNCombined.Final$Treatment,LNCombined.Final$cellType))
colnames(MajorTypeStatus)=c("Treatment","CellType","Number")
MajorTypeStatus$CellType=factor(MajorTypeStatus$CellType,levels=MajorTypeList)
g=ggplot(MajorTypeStatus,aes(x=Treatment, y=Number, fill=CellType)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=MajorTypeColors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("CellTypeIdentification/LNCombined.Final.Treatment.CellRatio.pdf",width=4,height=4)
print(g)
dev.off()


##############################################################################
############## Output cell type specific markers for Cibersoftx ##############
##############################################################################
setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/")
LNCombined.Final=readRDS("LN.Combined.Final.rds")
#32270 features across 463054 samples within 1 assay

CellTypeList=as.character(unique(LNCombined.Final$cellType))
CellType.list <- SplitObject(LNCombined.Final, split.by = "cellType")
CellType.list <- lapply(CellTypeList, FUN = function(x) {
    SubType=subset(LNCombined.Final,cellType==x)
    if(ncol(SubType)>5000){
       SubTypeDownSampling=SubType[, sample(colnames(SubType), size = 5000, replace=F)]
    }else{
        SubTypeDownSampling=SubType
    }
    return(SubTypeDownSampling)
 }
)

LNCombinedDownsampling=merge(CellType.list[[1]],CellType.list[-1])
saveRDS(LNCombinedDownsampling,file="LNCombinedDownsampling.rds")

CellType.markers <- FindAllMarkers(object = LNCombinedDownsampling,only.pos = TRUE,test.use = "wilcox",min.pct=0.25,logfc.threshold = 0.25)
write.table(CellType.markers,file="CellTypeIdentification/LN.CellType.Marker.txt",sep="\t",quote=F)


count=as.matrix(LNCombinedDownsampling@assays$RNA@counts)
all(colnames(count)==rownames(LNCombinedDownsampling@meta.data))
colnames(count)=LNCombinedDownsampling$cellType
GeneNames=data.frame("GeneSymbol"=rownames(count))
count.df=cbind(GeneNames,count)
MarkerGene.count=count.df[unique(CellType.markers$gene),]
write.table(MarkerGene.count,file="Cibersoftx/LN_SignatureMatrixFile.txt",row.names=F,sep="\t",quote=F)


############## Output cell type specific markers expression ##############
CellTypeExpr=AverageExpression(object = LNCombinedDownsampling)$RNA
CellTypeMarkers=read.table("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/CellTypeIdentification/LN.CellType.Marker.txt",header=T,row.names=1,sep="\t")
CellTypeMarkers=CellTypeMarkers[CellTypeMarkers$avg_log2FC>1,]
CellTypeMarkerExpr=CellTypeExpr[unique(CellTypeMarkers$gene),]

pdf("CellTypeIdentification/CellTypeMarker.heatmap.pdf",width=8,height=3)
pheatmap(t(CellTypeMarkerExpr),scale="column",clustering_method="ward.D2",show_colnames=F)
dev.off()

t=pheatmap(t(CellTypeMarkerExpr),scale="column",clustering_method="ward.D2",show_colnames=F)
dev.off()
GeneInfo=t$tree_col
GeneOrder=GeneInfo$label[GeneInfo$order]
CellTypeList=c("Macrophage","Epithelial cell","MKI67 cell","Germinal center B cell","B cell","Endothelial cell","Fibroblast cell","Neutrophial","NK cell","CD4 T cell","CD8 T cell","Exhausted T cell","Plasma cell","Mast cell")
pdf("CellTypeIdentification/CellTypeMarker.heatmap.modifed.pdf",width=8,height=4)
pheatmap(t(CellTypeMarkerExpr[GeneOrder,CellTypeList]),color=viridis(10),cluster_row=F,cluster_col=F,scale="column",clustering_method="ward.D2",show_colnames=F)
dev.off()


##############################################################################
############## subtype identification using cell type ########################
##############################################################################
library(pheatmap)
setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/")
LNCombined.Final=readRDS("LN.Combined.Final.rds")

SampleType=unique(LNCombined.Final@meta.data[,c("Sample","platform","Treatment")])
SampleType$LNStatus=ifelse(SampleType$Sample%in%c("CRC1_MLN1","CRC1_MLN2"),"Pos",ifelse(SampleType$Sample%in%c("CRC1_T","CRC2_T","CRC3_T"),"Tissue","Neg"))
SampleType$Batch=ifelse(SampleType$platform=="10XGenomics","Batch1",ifelse(SampleType$platform=="NuoHe","Batch2","Batch3"))
rownames(SampleType)=SampleType$Sample
SampleType$Sample=NULL
SampleType$platform=NULL

SampleTypeNumber=as.data.frame.matrix(table(LNCombined.Final$Sample,LNCombined.Final$cellType))
t=pheatmap(t(SampleTypeNumber),annotation_col=SampleType,clustering_method="ward.D2",scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
SampleAnno=data.frame(cluster=cutree(t$tree_col,k=5))
SampleType=SampleType[rownames(SampleAnno),]
SampleType$GroupByCellType=ifelse(SampleAnno$cluster%in%c(1,2),"NLN_C1",ifelse(SampleAnno$cluster%in%c(3,4),"PLN_C1","NLN_C2"))
SampleType$Treatment=ifelse(SampleType$Treatment%in%c("ICI+","ICI-"),SampleType$Treatment,"Not Available")
ann_colors = list(
                    GroupByCellType = c("PLN_C1"="Gold","NLN_C1"="Wheat","NLN_C2"="Orange"), 
                    Batch = c("Batch1" = "LightCyan", Batch2 = "Azure",Batch3="Thistle"),
                    LNStatus=c("Pos"="Orange","Neg"="LightCyan","Tissue"="LightGrey"),
                    Treatment=c("ICI+"="firebrick3","ICI-"="Bisque","Not Available"="LightGrey")
                )
pdf("LNSubTypeIdentification/LNSubTypeIdentification.heatmap.pdf",width=15,height=5)
pheatmap(t(SampleTypeNumber),annotation_col=SampleType[,c("GroupByCellType","LNStatus","Treatment")],annotation_colors=ann_colors,clustering_method="ward.D2",scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

SampleType$Sample=rownames(SampleType)
LNINfo=LNCombined.Final@meta.data
LNINfo$CellID=rownames(LNINfo)
LNINfo.merge=merge(LNINfo,SampleType,by="Sample")
rownames(LNINfo.merge)=LNINfo.merge$CellID
LNINfo.merge=LNINfo.merge[rownames(LNINfo),]
all(rownames(LNINfo.merge)==colnames(LNCombined.Final))
LNCombined.Final$GroupByCellType=LNINfo.merge$GroupByCellType
saveRDS(LNCombined.Final,file="LN.Combined.Final.rds")


MajorTypeList=c("CD4 T cell","CD8 T cell","NK cell","Exhausted T cell","B cell","Germinal center B cell","Plasma cell","MKI67 cell","Mast cell","Neutrophial","Macrophage","Epithelial cell","Fibroblast cell","Endothelial cell")
MajorTypeColors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(LNCombined.Final$cellType)))
names(MajorTypeColors)=MajorTypeList


t=pheatmap(t(SampleTypeNumber),clustering_method="ward.D2",scale="row")
LNOrder=t$tree_col$label[t$tree_col$order]

CellTypeNumberInSample=data.frame(table(LNCombined.Final$Sample,LNCombined.Final$cellType))
TotalNumberInSample=data.frame(table(LNCombined.Final$Sample))
CellTypeNumberInSample.ratio=merge(CellTypeNumberInSample,TotalNumberInSample,by="Var1")
colnames(CellTypeNumberInSample.ratio)=c("Sample","CellType","CellTypeNumber","TotalNumber")
CellTypeNumberInSample.ratio$Ratio=CellTypeNumberInSample.ratio$CellTypeNumber/CellTypeNumberInSample.ratio$TotalNumber
CellTypeNumberInSample.ratio$Sample=factor(CellTypeNumberInSample.ratio$Sample,levels=LNOrder)
CellTypeNumberInSample.ratio$CellType=factor(CellTypeNumberInSample.ratio$CellType,levels=MajorTypeList)
g=ggplot(CellTypeNumberInSample.ratio[CellTypeNumberInSample.ratio$CellType%in%c("B cell","Germinal center B cell","Fibroblast cell","Endothelial cell"),],aes(x=Sample, y=Ratio, fill=CellType)) +
      geom_bar(stat = "identity") +
      facet_grid(CellType~.,scale="free_y")+
      scale_fill_manual(values=MajorTypeColors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("LNSubTypeIdentification/LNCombined.Final.Treatment.CellRatio.BarPlot.pdf",width=14,height=6)
print(g)
dev.off()

g=ggplot(CellTypeNumberInSample.ratio,aes(x=Sample, y=Ratio, fill=CellType)) +
      geom_bar(stat = "identity") +
      facet_grid(CellType~.,scale="free_y")+
      scale_fill_manual(values=MajorTypeColors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("LNSubTypeIdentification/LNCombined.Final.Treatment.CellRatio.AllType.BarPlot.pdf",width=14,height=8)
print(g)
dev.off()


CellTypeNumberInSample.ratio.df=reshape2::dcast(CellTypeNumberInSample.ratio,CellType~Sample,value.var ="Ratio")
rownames(CellTypeNumberInSample.ratio.df)=CellTypeNumberInSample.ratio.df$CellType
CellTypeNumberInSample.ratio.df$CellType=NULL



sampleInfo=unique(LNCombined.Final@meta.data[,c("Sample","Type","Treatment","platform")])
CellInfo=LNCombined.Final@meta.data
t=sapply(unique(CellInfo$Sample), function(s) mean(as.numeric(CellInfo[CellInfo$Sample==s,"nFeature_RNA"])))
all(names(t)==sampleInfo$Sample)
sampleInfo$nFeature_RNA=t
t=sapply(unique(CellInfo$Sample), function(s) mean(as.numeric(CellInfo[CellInfo$Sample==s,"nCount_RNA"])))
all(names(t)==sampleInfo$Sample)
sampleInfo$nCount_RNA=t
write.table(sampleInfo,file="SampleInfo.txt",sep="\t",quote=F,row.names=F)



##############################################################################
############## End and Fibroblast ########################
##############################################################################
setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/")
LNCombined.Final=readRDS("LN.Combined.Final.rds")
EndAndFib=subset(LNCombined.Final,cellType%in%c("Fibroblast cell","Endothelial cell"))

EndAndFib <- NormalizeData(EndAndFib, normalization.method = "LogNormalize", scale.factor = 10000)
EndAndFib <- FindVariableFeatures(EndAndFib, selection.method = "vst", nfeatures = 2000)
EndAndFib <- ScaleData(EndAndFib)
EndAndFib <- RunPCA(EndAndFib, features = VariableFeatures(object = EndAndFib))
EndAndFib <- RunHarmony(EndAndFib, group.by.vars=c("platform"))
EndAndFib<- RunUMAP(EndAndFib, reduction = "harmony", dims = 1:10)
EndAndFib<- FindNeighbors(EndAndFib, reduction = "harmony", dims = 1:10)
EndAndFib<- FindClusters(EndAndFib, resolution = 2)

pdf("EndAndFib/EndAndFib.Cluster.pdf",width=6,height=5)
DimPlot(EndAndFib,label=T,group.by="seurat_clusters")&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("EndAndFib/EndAndFib.CellType.pdf",width=7,height=5)
DimPlot(EndAndFib,label=T,group.by="cellType")&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("EndAndFib/EndAndFib.CellMarker.pdf",width=10,height=10)
FeaturePlot(EndAndFib,pt.size=3,raster=TRUE,features = c("LYZ","PROX1","RBP7","KDR","ACKR1","FCAMR","CR2","CR1","CRISPLD2","NTRK3","PLN","ACTA2","ADH1B","C3","COL5A1","VCAN","PDPN","IL7R"))&NoLegend()&NoAxes()&theme(plot.title=element_text(face="italic"),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("EndAndFib/EndAndFib.CellMarker.AllCell.pdf",width=10,height=10)
FeaturePlot(LNCombined.Final,raster=TRUE,features = c("COL1A1","ENG","LYZ","PROX1","CA4","KDR","ACKR1","CXCL13","CRISPLD2","NTRK3","PLN","HIGD1B","ACTA2","ADH1B","PDPN","IL7R"))&NoLegend()&NoAxes()&theme(plot.title=element_text(face="italic"),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("EndAndFib/EndAndFib.CellMarker.pdf",width=10,height=12)
Plot_Density_Custom(seurat_object =EndAndFib,features = rev(c("LYZ","PROX1","RBP7","KDR","ACKR1","FCAMR","CR2","CR1","CRISPLD2","NTRK3","PLN","ACTA2","ADH1B","C3","COL5A1","VCAN","PDPN","IL7R")))&NoLegend()&NoAxes()&theme(plot.title=element_text(face="italic"),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()



EndAndFib.markers <- FindAllMarkers(EndAndFib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(EndAndFib.markers,file="EndAndFib/Subcluster.Marker.txt",sep="\t",quote=F)

EndAndFib.markers=read.table("EndAndFib/Subcluster.Marker.txt",sep="\t",header=T)
for(i in unique(EndAndFib$seurat_clusters)){
t=FeaturePlot(EndAndFib, features = EndAndFib.markers[EndAndFib.markers$cluster==i,"gene"][1:16])&NoLegend()&NoAxes()&theme(plot.title=element_text(face="italic"),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
tiff(paste0("EndAndFib/SubSclusterMarker/C",i,"_Marker_Harmony.Subcluster.tiff",sep=""),width=1000,height=1000)
print(t)
dev.off()
}

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
pdf("EndAndFib/Harmony.nFeaturesRNA.pdf",height=8,width=15)
VlnPlot(EndAndFib,group.by = "seurat_clusters", features =feats, pt.size = 0, ncol = 1) +NoLegend()
dev.off()

Idents(EndAndFib)=EndAndFib$seurat_clusters
new.cluster.ids<-c("Fib IL7R","Fib IL7R","aBEC","cBEC","FRC","SMC","FRC","Fib IL7R","SMC","vBEC","fLEC","SMC","Fib ADH1B","Fib ADH1B","vBEC","Fib NTRK3","Fib ADH1B","LEC","SMC","FDC","FDC","cBEC","aBEC","FRC","aBEC")
names(new.cluster.ids) <- levels(EndAndFib)
EndAndFib <- RenameIdents(EndAndFib, new.cluster.ids)
EndAndFib$SubCellType=Idents(EndAndFib)
EndAndFib$SubCellType=factor(EndAndFib$SubCellType,levels=c("Fib IL7R","FRC","Fib ADH1B","SMC","Fib NTRK3","FDC","vBEC","aBEC","cBEC","LEC","fLEC"))
SubCellTypeColors <- c(colorRampPalette(brewer.pal(6, "Dark2"))(6),colorRampPalette(brewer.pal(6, "Set3"))(5))
pdf("EndAndFib/EndAndFib.SubCellType.pdf",width=5,height=3.5)
DimPlot(EndAndFib,cols =SubCellTypeColors,group.by="SubCellType")&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()

EndAndFib.CellType.markers <- FindAllMarkers(EndAndFib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(EndAndFib.CellType.markers,file="EndAndFib/EndAndFib.CellType.markers",sep="\t",quote=F)

######Clustered_DotPlot ###########
all_markers <- EndAndFib.CellType.markers %>% Add_Pct_Diff()
top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7, named_vector = FALSE, make_unique = TRUE)
pdf("EndAndFib/EndAndFib.SubCellType.markers.Clustered_DotPlot.pdf",width=8,height=8)
Clustered_DotPlot(seurat_object = EndAndFib, features = top_markers)
dev.off()

######Stacked_VlnPlot ###########
MarkerGene=rev(c("LYZ","PROX1","RBP7","KDR","ACKR1","FCAMR","CR2","CR1","CRISPLD2","NTRK3","PLN","ACTA2","ADH1B","C3","COL5A1","VCAN","PDPN","IL7R"))
pdf("EndAndFib/EndAndFib.Final.CellType.StackedVlnPlot.pdf",height=8,width=6)
Stacked_VlnPlot(seurat_object = EndAndFib,plot_spacing =0,group.by="SubCellType", features = MarkerGene, x_lab_rotate = 90,colors_use = SubCellTypeColors)
dev.off()

###### DotPlot in All cells###########
LNCombined.Final$cellType=factor(LNCombined.Final$cellType,levels=c("CD4 T cell","CD8 T cell","NK cell","Exhausted T cell","B cell","Germinal center B cell","Plasma cell","MKI67 cell","Mast cell","Neutrophial","Macrophage","Epithelial cell","Fibroblast cell","Endothelial cell"))
CellTypeColors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(LNCombined.Final$cellType)))

MarkerGene=rev(c("FCAMR","CR2","CR1","KDR","COL5A1","VCAN","PDPN"))
pdf("EndAndFib/EndAndFib.Final.CellType.Dotplot.InAllType.pdf",width=6,height=5)
DotPlot(LNCombined.Final,features =MarkerGene,group.by="cellType")&RotatedAxis()&theme(axis.title.y=element_blank(),axis.title.x=element_blank(),plot.title=element_text(face="italic"),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()

saveRDS(EndAndFib,file="EndAndFib/EndAndFib.Final.rds")


##############################################################################
############## cell type variation between pos and neg ########################
##############################################################################
setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/")
EndAndFib=readRDS("EndAndFib/EndAndFib.Final.rds")

SubCellTypeColors <- c(colorRampPalette(brewer.pal(6, "Dark2"))(6),colorRampPalette(brewer.pal(6, "Set3"))(5))
pdf("LNSubTypeVariation/EndAndFib.Variation.pdf",width=5,height=3.5)
DimPlot(EndAndFib,cols =SubCellTypeColors,group.by="SubCellType",split.by="GroupByCellType",ncol=2)&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()


EndAndFib$GroupByCellType=ifelse(EndAndFib$Sample%in%c("CRC1_T","CRC2_T","CRC3_T"),"Tissue",EndAndFib$GroupByCellType)
CellTypeGroup=data.frame(table(EndAndFib$GroupByCellType,EndAndFib$SubCellType))
colnames(CellTypeGroup)=c("Group","SubCellType","Number")

FibGroup=CellTypeGroup[CellTypeGroup$SubCellType%in%c("Fib IL7R","FRC","Fib C3","SMC","Fib NTRK3","FDC"),]
FibCellType_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(FibGroup$SubCellType)))
FibGroup$SubCellType=factor(FibGroup$SubCellType,levels=c(sort(unique(FibGroup$SubCellType))))
g=ggplot(FibGroup,aes(x=Group, y=Number, fill=SubCellType)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=FibCellType_colors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("LNSubTypeVariation/FibSubTypeVariation.pdf",width=4,height=4)
print(g)
dev.off()

EndGroup=CellTypeGroup[CellTypeGroup$SubCellType%in%c("vBEC","aBEC","cBEC","LEC","fLEC"),]
EndCellType_colors <- colorRampPalette(brewer.pal(8, "Set3"))(length(unique(EndGroup$SubCellType)))
EndGroup$SubCellType=factor(EndGroup$SubCellType,levels=c(sort(unique(EndGroup$SubCellType))))
g=ggplot(EndGroup,aes(x=Group, y=Number, fill=SubCellType)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=EndCellType_colors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("LNSubTypeVariation/EndSubTypeVariation.pdf",width=4,height=4)
print(g)
dev.off()


##############################################################################
############## cell type variation between Treatment  ########################
##############################################################################

EndAndFib.Treatment=subset(EndAndFib,Treatment%in%c("ICI+","ICI-"))
CellTypeGroup=data.frame(table(EndAndFib.Treatment$Treatment,EndAndFib.Treatment$SubCellType))
colnames(CellTypeGroup)=c("Treatment","SubCellType","Number")

SubCellTypeColors <- c(colorRampPalette(brewer.pal(6, "Dark2"))(6),colorRampPalette(brewer.pal(6, "Set3"))(5))
pdf("LNSubTypeVariation/EndAndFib.Treatment.Variation.pdf",width=7,height=3)
DimPlot(EndAndFib.Treatment,cols =SubCellTypeColors,group.by="SubCellType",split.by="Treatment",ncol=2)&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()

FibGroup=CellTypeGroup[CellTypeGroup$SubCellType%in%c("Fib IL7R","FRC","Fib C3","SMC","Fib NTRK3","FDC"),]
FibCellType_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(FibGroup$SubCellType)))
FibGroup$SubCellType=factor(FibGroup$SubCellType,levels=c(sort(unique(FibGroup$SubCellType))))
FibGroup$Treatment=factor(FibGroup$Treatment,levels=c("ICI+","ICI-"))
g=ggplot(FibGroup,aes(x=Treatment, y=Number, fill=SubCellType)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=FibCellType_colors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("LNSubTypeVariation/Treatment_FibSubTypeVariation.pdf",width=3,height=3)
print(g)
dev.off()


EndGroup=CellTypeGroup[CellTypeGroup$SubCellType%in%c("vBEC","aBEC","cBEC","LEC","fLEC"),]
EndCellType_colors <- colorRampPalette(brewer.pal(8, "Set3"))(length(unique(EndGroup$SubCellType)))
EndGroup$SubCellType=factor(EndGroup$SubCellType,levels=c(sort(unique(EndGroup$SubCellType))))
EndGroup$Treatment=factor(EndGroup$Treatment,levels=c("ICI+","ICI-"))
g=ggplot(EndGroup,aes(x=Treatment, y=Number, fill=SubCellType)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=EndCellType_colors) +
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("LNSubTypeVariation/Treatment_EndSubTypeVariation.pdf",width=3,height=3)
print(g)
dev.off()


##############################################################################
############## cell type variation (DEG) between Treatment  ##################
##############################################################################
library(pheatmap)
setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/")
EndAndFib=readRDS("EndAndFib/EndAndFib.Final.rds")
EndAndFib.Treatment=subset(EndAndFib,Treatment%in%c("ICI+","ICI-"))
Idents(EndAndFib.Treatment)=EndAndFib.Treatment$Treatment

Fib.Treatment=subset(EndAndFib.Treatment,SubCellType%in%c("Fib IL7R","FRC","Fib C3","SMC","Fib NTRK3","FDC"))
FibDEG=FindMarkers(Fib.Treatment,ident.1="ICI+",ident.2="ICI-",test.use = "MAST",latent.vars ="platform",logfc.threshold = 0.25)
SigDEG=FibDEG[FibDEG$p_val<0.01,]
SigDEG$Pattern=ifelse(SigDEG$avg_log2FC>0,"Up","Down")
table(SigDEG$Pattern)
write.table(FibDEG,file="LNSubTypeVariation/Treatment_EndAndFib_DEG.txt",sep="\t",quote=F)

GeneAnno=SigDEG[,c("Pattern","p_val")]
GeneAnno$logP=-log10(GeneAnno$p_val)
GeneAnno$p_val=NULL
SubTypeExpr=AverageExpression(object =Fib.Treatment,group.by="SubCellType")$RNA
SigGenes=intersect(rownames(SigDEG),rownames(SubTypeExpr))
length(SigGenes)


t=pheatmap(SubTypeExpr[rownames(SigDEG),],scale="row",clustering_method="average")
dev.off()
geneInfo=data.frame(Cluster=cutree(t$tree_row,k=6))
all(rownames(geneInfo)==rownames(GeneAnno))
GeneAnno$G=paste("C",geneInfo$Cluster,sep="")
GeneAnno[GeneAnno$G%in%c("C1"),"Group"]="Fib NTRK3"
GeneAnno[GeneAnno$G%in%c("C4","C6"),"Group"]="FDC"
GeneAnno[GeneAnno$G%in%c("C5"),"Group"]="Fib IL7R"
GeneAnno[GeneAnno$G%in%c("C2"),"Group"]="FRC"
GeneAnno[GeneAnno$G%in%c("C3"),"Group"]="SMC"
GeneAnno$G=NULL
GeneOrder=t$tree_row$label[t$tree_row$order]
CellTypeList=c("Fib NTRK3","FDC","Fib IL7R","FRC","SMC")
GroupColors=colorRampPalette(brewer.pal(6, "Set3"))(length(unique(GeneAnno$Group)))
names(GroupColors)=unique(GeneAnno$Group)
ann_colors = list(
                    Pattern = c("Up"="Red","Down"="Blue"), 
                    logP = c("white","Purple"),
                    Group=GroupColors
                )
pdf("LNSubTypeVariation/Treatment_Fib_DEG_CellTypeExpr_Heatmap.pdf",width=5,height=5)
pheatmap(SubTypeExpr[GeneOrder,CellTypeList],cluster_row=F,cluster_col=F,annotation_colors=ann_colors,scale="row",show_rownames=F,annotation_row=GeneAnno,clustering_method="average", color = colorRampPalette(c("LightBlue", "white", "Orange"))(50))
dev.off()

GeneAnnoNumber=data.frame(table(GeneAnno$Group,GeneAnno$Pattern))
colnames(GeneAnnoNumber)=c("CellType","Pattern","Number")
g=ggplot(GeneAnnoNumber,aes(x=CellType, y=Number, fill=Pattern)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values=c("Blue","Red")) +
      theme_bw()+
      coord_flip()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("LNSubTypeVariation/Treatment_FibEndAllDEG_Number_Barplot.pdf",width=4,height=2)
print(g)
dev.off()

write.table(GeneAnno[order(GeneAnno$Group),],file="LNSubTypeVariation/Treatment_Fib_DEG_CellType.txt",sep="\t",quote=F)


targetGenes=c("PDPN","COL1A1","VCAN","FDCSP","FCAMR","CR1","CR2")
t=DotPlot(subset(EndAndFib.Treatment,SubCellType=="FDC"),group.by="Treatment",features=targetGenes)&RotatedAxis()
graphData=t$data
graphData$id=factor(graphData$id,levels=c("ICI+","ICI-"))
g=ggplot(graphData,aes(x=id,y=features.plot,color=avg.exp.scaled,size=pct.exp))+ 
geom_point()+
theme_bw()+
scale_colour_gradient(low="blue",high="red")+
theme(axis.title.x=element_blank(),axis.title.y=element_blank())
pdf("LNSubTypeVariation/Treatment_FDCAndFDC_MarkgeGenesVariation_Point_FDCOnly.pdf",width=3,height=2)
print(g)
dev.off()

targetGenes=c("SPIB","LTF","RELB")
t=DotPlot(EndAndFib.Treatment,group.by="Treatment",features=rev(targetGenes))&RotatedAxis()
graphData=t$data
graphData$id=factor(graphData$id,levels=c("ICI+","ICI-"))
g=ggplot(graphData,aes(x=id,y=features.plot,color=avg.exp.scaled,size=pct.exp))+ 
geom_point()+
theme_bw()+
scale_colour_gradient(low="blue",high="red")+
theme(axis.title.x=element_blank(),axis.title.y=element_blank())
pdf("LNSubTypeVariation/Treatment_FDC_TFVariation_Point.pdf",width=3,height=1.5)
print(g)
dev.off()

targetGenes=read.table("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/RELB.targetGene.txt")
targetGenes=targetGenes[,1]
t=DotPlot(EndAndFib.Treatment,group.by="Treatment",features=rev(targetGenes))&RotatedAxis()
graphData=t$data
graphData$id=factor(graphData$id,levels=c("ICI+","ICI-"))
g=ggplot(graphData,aes(x=id,y=features.plot,color=avg.exp.scaled,size=pct.exp))+ 
geom_point()+
theme_bw()+
scale_colour_gradient(low="blue",high="red")+
theme(axis.title.x=element_blank(),axis.title.y=element_blank())
pdf("LNSubTypeVariation/Treatment_RELBTArgetGenes_Point.pdf",width=3,height=10)
print(g)
dev.off()

targetGenes=read.table("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/SPIB.targetGene.txt")
targetGenes=targetGenes[,1]
t=DotPlot(EndAndFib.Treatment,group.by="Treatment",features=rev(targetGenes))&RotatedAxis()
graphData=t$data
graphData$id=factor(graphData$id,levels=c("ICI+","ICI-"))
g=ggplot(graphData,aes(x=id,y=features.plot,color=avg.exp.scaled,size=pct.exp))+ 
geom_point()+
theme_bw()+
scale_colour_gradient(low="blue",high="red")+
theme(axis.title.x=element_blank(),axis.title.y=element_blank())
pdf("LNSubTypeVariation/Treatment_SPIBTArgetGenes_Point.pdf",width=3,height=10)
print(g)
dev.off()


#########################################################################################
#########################################################################################
EndAndFib.Treatment=subset(EndAndFib,Treatment%in%c("ICI+","ICI-"))
Idents(EndAndFib.Treatment)=EndAndFib.Treatment$Treatment
table(EndAndFib.Treatment$platform,EndAndFib.Treatment$Sample)
#check the sample in each platform
#EndAndFib.Treatment.tmp=subset(EndAndFib.Treatment,platform%in%c("NuoHe"))
EndAndFib.Treatment.tmp=EndAndFib.Treatment
#remove dMMR sample
pMMRSample=c("C4_LN1","C4_LN2","C4_LN3","C6_LN1","C6_LN2","C6_LN3","C6_LN4")
dMMRSample=setdiff(unique(EndAndFib.Treatment.tmp$Sample),pMMRSample)
EndAndFib.Treatment.tmp=subset(EndAndFib.Treatment.tmp,Sample%in%dMMRSample)
unique(EndAndFib.Treatment.tmp@meta.data[,c("platform","Sample","Treatment")])
#                          platform Sample Treatment
#C1_LN1_GATCGCGAGACAATAC-1    NuoHe C1_LN1      ICI+
#C1_LN2_ACGGGTCAGAATGTGT-1    NuoHe C1_LN2      ICI+
#C1_LN3_AGTGTCATCGATAGAA-1    NuoHe C1_LN3      ICI+
#C2_LN1_AACTCCCGTCGCATCG-1    NuoHe C2_LN1      ICI+
#C2_LN3_AAATGCCCACGGACAA-1    NuoHe C2_LN3      ICI+
#C3_LN1_AAGGCAGTCCGTCAAA-1    NuoHe C3_LN1      ICI+
#C3_LN2_ACACCGGTCCAGTATG-1    NuoHe C3_LN2      ICI+
#C3_LN3_AAAGTAGAGCTAGTTC-1    NuoHe C3_LN3      ICI+
#C5_LN1_TCGTAGAAGGGATCTG-1    NuoHe C5_LN1      ICI-
#C5_LN2_TGGTTCCAGTCCTCCT-1    NuoHe C5_LN2      ICI-
#C5_LN3_ACGGGCTCACGTAAGG-1    NuoHe C5_LN3      ICI-
#C5_LN4_CTAGTGAAGAAGGCCT-1    NuoHe C5_LN4      ICI-
#C7_LN1_CACCAGGAGAAGGGTA-1    NuoHe C7_LN1      ICI-
#C7_LN2_ACGATGTTCCTACAGA-1    NuoHe C7_LN2      ICI-
#C7_LN3_AACTGGTCATGCTAGT-1    NuoHe C7_LN3      ICI-
#C8_LN1_AAGCCGCCACAGGAGT-1    NuoHe C8_LN1      ICI+
#C8_LN2_ACATACGGTCCAAGTT-1    NuoHe C8_LN2      ICI+
#C8_LN3_CAGTCCTTCATGGTCA-1    NuoHe C8_LN3      ICI+
#C8_LN4_ACGCAGCTCCGATATG-1    NuoHe C8_LN4      ICI+
#C9_LN1_AATCCAGAGTTTAGGA-1    NuoHe C9_LN1      ICI+
#C9_LN2_AAGGAGCTCTTCGAGA-1    NuoHe C9_LN2      ICI+
#C9_LN3_AACTCAGCACATGTGT-1    NuoHe C9_LN3      ICI+
#C9_LN4_AAGGTTCCATGGAATA-1    NuoHe C9_LN4      ICI+

targetGenes=c("ACTA2","PDPN","COL1A1","VCAN","FDCSP","FCAMR","CR1","CR2")
t=DotPlot(EndAndFib.Treatment.tmp,group.by="Treatment",features=targetGenes)&RotatedAxis()
graphData=t$data
graphData$id=factor(graphData$id,levels=c("ICI+","ICI-"))
g=ggplot(graphData,aes(x=id,y=features.plot,color=avg.exp.scaled,size=pct.exp))+ 
geom_point()+
theme_bw()+
scale_colour_gradient(low="blue",high="red")+
theme(axis.title.x=element_blank(),axis.title.y=element_blank())
pdf("LNSubTypeVariation/Treatment_FDCAndFDC_MarkgeGenesVariation_Point_FibAndEnd_OnlyIn_pMMRSample.pdf",width=3,height=2.5)
print(g)
dev.off()

#########################################################################################
#########################################################################################


FibDEG=FindMarkers(Fib.Treatment,ident.1="ICI+",ident.2="ICI-",test.use = "MAST",latent.vars ="platform",logfc.threshold = 0)

RELB.targetGen=read.table("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/RELB.targetGene.txt")
RELB.targetGen=RELB.targetGen[,1]
RELB.targetGen.DEG=FibDEG[intersect(RELB.targetGen,rownames(FibDEG)),]
library(ggrepel)
hs_data=data.frame(FibDEG)
hs_data$threshold = as.factor(ifelse(hs_data$p_val>0.05,"NoSig",ifelse(abs(hs_data$avg_log2FC)>0.25,ifelse(hs_data$avg_log2FC>0.25,"SigUp","SigDown"),ifelse(hs_data$avg_log2FC>0,"Up","Down"))))
table(hs_data$threshold)
#Down   NoSig SigDown   SigUp      Up 
#4420   21760    1871    2959    3539
hs_data$threshold=factor(hs_data$threshold,levels=c("SigDown","Down","NoSig","Up","SigUp"))
hs_data$FibDEG = ifelse(rownames(hs_data)%in%RELB.targetGen,"RELB Target genes","Others")
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = avg_log2FC, y = -log10(p_val), size=log2(pct.1+1),colour=FibDEG, label =ID )) +
  geom_point(alpha=0.9) +
  theme_bw() + 
  scale_size(range = c(1, 1))+
  scale_color_manual(values=c("LightGrey","red")) +
  xlim(c(-2, 2)) + #ylim(0,25)+
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
    geom_text_repel(
    data =hs_data[intersect(RELB.targetGen,rownames(hs_data)),],
    aes(label = ID),
    size = 3,
    max.overlaps=5,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/RELB.targetGene_Volcano_NoLabel.pdf",width=5,height=4)
print(t)
dev.off()



library(presto) #download zip file from github, unzip, install.packages(path,rep=NULL,type="resource")
library(dplyr)
library(msigdbr)
library(tibble)
library(fgsea)
setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/")
EndAndFib=readRDS("EndAndFib/EndAndFib.Final.rds")
EndAndFib.Treatment=subset(EndAndFib,Treatment%in%c("ICI+","ICI-"))
Idents(EndAndFib.Treatment)=EndAndFib.Treatment$Treatment
Fib.Treatment=subset(EndAndFib.Treatment,SubCellType%in%c("Fib IL7R","FRC","Fib C3","SMC","Fib NTRK3","FDC"))
Treatment.DEG <- wilcoxauc(Fib.Treatment, 'Treatment')
table(Treatment.DEG$group)

RELB.targetGene=read.table("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/RELB.targetGene.txt")
SPIB.targetGene=read.table("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/SPIB.targetGene.txt")
LTF.targetGene=read.table("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/LTF.targetGene.txt")
fgsea_sets=list(RELBTarget=RELB.targetGene[,1],SPIBTarget=SPIB.targetGene[,1],LTFTarget=LTF.targetGene[,1])
clusterCell<- Treatment.DEG %>% dplyr::filter(group == "ICI+") %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
ranks=na.omit(ranks)
fgseaRes <- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0, nPermSimple = 10000)
fwrite(fgseaRes, file=paste0("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/TargetGeneBetweenTreatment.txt",sep=""), sep="\t", sep2=c("", " ", ""),quote=F)

pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/ICITreatment_RELBTarget.pdf",width=3,height=2)
plotEnrichment(fgsea_sets[["RELBTarget"]],ranks)
dev.off()
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/ICITreatment_SPIBTarget.pdf",width=3,height=2)
plotEnrichment(fgsea_sets[["SPIBTarget"]],ranks)
dev.off()


CellTypeGroup=data.frame(table(EndAndFib.Treatment$Treatment,EndAndFib.Treatment$SubCellType))
colnames(CellTypeGroup)=c("Treatment","SubCellType","Number")
LowCellNumberType=CellTypeGroup[CellTypeGroup$Number<50,]

TargetCellType=setdiff(unique(CellTypeGroup$SubCellType),unique(LowCellNumberType$SubCellType))
for(type in TargetCellType){
    SubType=subset(EndAndFib.Treatment,SubCellType==type)
    if(type==TargetCellType[1]){
        AllDEG=FindMarkers(SubType,ident.1="ICI+",ident.2="ICI-",test.use = "MAST",latent.vars ="platform",logfc.threshold = 0.25)
        AllDEG$CellType=type
        AllDEG$Gene=rownames(AllDEG)
    }else{
        DEG=FindMarkers(SubType,ident.1="ICI+",ident.2="ICI-",test.use = "MAST",latent.vars ="platform",logfc.threshold = 0.25)
        DEG$CellType=type
        DEG$Gene=rownames(DEG)
        AllDEG=rbind(AllDEG,DEG)
    }
}
write.table(AllDEG,file="LNSubTypeVariation/Treatment_SubType_DEG.txt",sep="\t",quote=F)
SigDEG=AllDEG[AllDEG$p_val<0.01,]
SigDEG$Pattern=ifelse(SigDEG$avg_log2FC>0,"Up","Down")
SigDEGNumber=data.frame(table(SigDEG$CellType,SigDEG$Pattern))
colnames(SigDEGNumber)=c("CellType","Pattern","Number")

g=ggplot(SigDEGNumber,aes(x=CellType, y=Number, fill=Pattern)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values=c("Blue","Red")) +
      theme_bw()+
      coord_flip()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("LNSubTypeVariation/Treatment_SubType_DEG_Number_Barplot.pdf",width=4,height=2)
print(g)
dev.off()




