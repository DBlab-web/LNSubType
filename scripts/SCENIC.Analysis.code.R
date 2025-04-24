library(SCENIC)
library(Seurat)
library(RColorBrewer)


setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/")
EndAndFib=readRDS("EndAndFib/EndAndFib.Final.rds")

library(data.table)
devtools::install_local("/groups/g900008/home/share/SCENIC-1.1.2.zip")
#https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/
#old version
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
scenicOptions <- initializeScenic(org="hgnc",dbDir="/lkn_lab/bixiaoman/Software/SCENIC/cisTarget_databases", nCores=10)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC")
### Co-expression network
exprMat <- as.matrix(EndAndFib@assays$RNA$data)
genesKept <- geneFiltering(exprMat, scenicOptions)
length(genesKept)
#9590
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
source("/groups/g900008/home/share/runSCENIC_2_createRegulons.R")
exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
#Number of links between TFs and targets (weight>=0.001): 3408251
#             [,1]
#nTFs         1014
#nTargets     9590
#nGeneSets    8097
#nLinks    4573458
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
#object 'indexCol' not found.#https://github.com/aertslab/SCENIC/issues/290
#Scoring database:  [Source file: hg19-500bp-upstream-7species.mc9nr.feather]
#Scoring database:  [Source file: hg19-tss-centered-10kb-7species.mc9nr.feather]
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)


CellTypeMarkers=read.table("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/EndAndFib.CellType.markers",header=T)

##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
#scRNAauc <- AddMetaData(scRNA, AUCmatrix)
#scRNAauc@assays$integrated <- NULL
#saveRDS(scRNAauc,'scRNAauc.rds')
cellInfo=EndAndFib@meta.data
all(rownames(AUCmatrix)==rownames(cellInfo))
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$SubCellType), function(cells) colMeans(AUCmatrix[cells,]))
dim(regulonActivity_byCellType)
TFTmp=grep("_extended",rownames(regulonActivity_byCellType),value=T)
regulonActivity_byCellType.filter=regulonActivity_byCellType[setdiff(rownames(regulonActivity_byCellType),TFTmp),]

TFList=data.frame(TFInfo=rownames(regulonActivity_byCellType.filter))
TFList$TFName=do.call(rbind,strsplit(as.character(TFList$TFInfo),split =" "))[,1]
TFList=TFList[TFList$TFName%in%CellTypeMarkers$gene,]
#59
regulonActivity_byCellType.filter=regulonActivity_byCellType.filter[TFList$TFInfo,]
regulonActivity_byCellType.filter=regulonActivity_byCellType.filter[Biobase::rowMax(regulonActivity_byCellType.filter)>0.01,]
dim(regulonActivity_byCellType.filter)

pdf("regulonActivity_byCellType.pdf",height=6,width=4)
pheatmap(regulonActivity_byCellType.filter,scale="row",show_rownames=F,clustering_method="ward.D2",color = viridis(10))
dev.off()
pdf("regulonActivity_byCellType_Name.pdf",height=10,width=4)
pheatmap(regulonActivity_byCellType.filter,scale="row",show_rownames=T,clustering_method="ward.D2",color = viridis(10))
dev.off()




setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/")
LNCombined.Final=readRDS("LN.Combined.Final.rds")
EndAndFib=readRDS("EndAndFib/EndAndFib.Final.rds")
all(colnames(EndAndFib)%in%colnames(LNCombined.Final))
#TRUE

table(EndAndFib$SubCellType)
table(LNCombined.Final$cellType)
LNAllInfo=LNCombined.Final@meta.data
FibEndInfo=EndAndFib@meta.data
LNAllInfo$CellID=rownames(LNAllInfo)
FibEndInfo$CellID=rownames(FibEndInfo)
LNAllInfo.merge=merge(LNAllInfo[,c("CellID","cellType")],FibEndInfo[,c("CellID","SubCellType")],by="CellID",all.x=TRUE)
rownames(LNAllInfo.merge)=LNAllInfo.merge$CellID
LNAllInfo.merge=LNAllInfo.merge[rownames(LNAllInfo),]
all(rownames(LNAllInfo)==rownames(LNAllInfo.merge))
LNAllInfo.merge$ModifiedSubCellType=ifelse(is.na(LNAllInfo.merge$SubCellType),as.character(LNAllInfo.merge$cellType),as.character(LNAllInfo.merge$SubCellType))
table(LNAllInfo.merge$ModifiedSubCellType)
LNCombined.Final$SubCellType=LNAllInfo.merge$ModifiedSubCellTyp
LNCombined.Final$SubCellType=factor(LNCombined.Final$SubCellType,levels=c("CD4 T cell","CD8 T cell","NK cell","Exhausted T cell","B cell","Germinal center B cell","Plasma cell","MKI67 cell","Mast cell","Neutrophial","Macrophage","Epithelial cell","Fib IL7R","FRC","Fib ADH1B","SMC","Fib NTRK3","FDC","vBEC","aBEC","cBEC","LEC","fLEC"))
table(LNCombined.Final$SubCellType)

targetTF=c("EGR1","CEBPD","SOX4","HES1","KLF10","HIF1A","ELF3","PRDM1","ETS2","XBP1","CEBPB","FOXF2","NFIC","FOXF1","PRNP","SMARCA1","NFIX","STAT1","SPIB","LTF","RELB","MEF2C","KLF4","ERG","SMAD1","MECOM","SOX18")
t=DotPlot(LNCombined.Final,group.by="SubCellType",features=rev(targetTF))&RotatedAxis()
graphData=t$data
graphData$CellType=ifelse(graphData$features.plot%in%c("STAT1","SPIB","LTF","RELB"),"FDC",ifelse(graphData$features.plot%in%c("MEF2C","KLF4","ERG","SMAD1","MECOM","SOX18"),"BEC","FRC"))
graphData$CellType=factor(graphData$CellType,levels=c("FRC","FDC","BEC"))
g=ggplot(graphData,aes(x=id,y=features.plot,color=avg.exp.scaled,size=pct.exp))+ 
geom_point()+
theme_bw()+
facet_grid(CellType~.,scale="free_y",space="free_y")+
scale_colour_gradient2(low="blue",mid="white",high="red",midpoint=0)+
theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"),axis.title.x=element_blank(),,axis.title.y=element_blank())
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/FRC_FDC_DotPlot.pdf",width=6,height=6)
print(g)
dev.off()


MarkerGene=c("CDKN1A","CDKN2A","CDKN2B")
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SenescenceMarker_VlnPlot.pdf",height=4,width=6)
Stacked_VlnPlot(seurat_object =LNCombined.Final,plot_spacing =0,group.by="SubCellType", features = MarkerGene, x_lab_rotate = 90)
dev.off()

pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SenescenceMarker_DotPlot.pdf",width=6,height=5)
DotPlot(LNCombined.Final,features =MarkerGene,group.by="SubCellType")&RotatedAxis()&theme(axis.title.y=element_blank(),axis.title.x=element_blank(),plot.title=element_text(face="italic"),panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()

#### Bulk RNA  ##################
setwd("D:/06_CRC/bulkRNAseq/BeiJingSample/")
LN_TPM_BeiJing=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/LN_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
LN_Info_BeiJing=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/LN_SubType_Info.txt",header=T,row.names=1,sep="\t")
LN_Info_BeiJing$GroupByGene=factor(LN_Info_BeiJing$GroupByGene,levels=sort(unique(LN_Info_BeiJing$GroupByGene)))
LN_Info_BeiJing=LN_Info_BeiJing[order(LN_Info_BeiJing$GroupByGene),]
LN_TPM_BeiJing=LN_TPM_BeiJing[,rownames(LN_Info_BeiJing)]
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info_BeiJing$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info_BeiJing$GroupByGene))

targetTF=c("EGR1","CEBPD","SOX4","HES1","KLF10","HIF1A","ELF3","PRDM1","ETS2","XBP1","CEBPB","FOXF2","NFIC","FOXF1","PRNP","SMARCA1","NFIX","STAT1","SPIB","LTF","RELB","MEF2C","KLF4","ERG","SMAD1","MECOM","SOX18")
targetTFExpr=LN_TPM_BeiJing[targetTF,]
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info_BeiJing$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info_BeiJing$GroupByGene))
ann_colors = list(
                    GroupByGene = GroupByGene_Colors, 
                    LNGroup=c("Pos"="Orange","Neg"="LightBlue")
                )
pdf("D:/06_CRC/Graph/Part4/TFExprInBeiJing.pdf",width=8,height=4)
pheatmap(targetTFExpr,annotation_col=LN_Info_BeiJing[,c("GroupByGene","LNGroup")],gaps_row = c(17,21),
    cluster_col=F,cluster_row=F,clustering_method="ward.D2",scale="row",annotation_colors = ann_colors,
    show_colnames=F,color = colorRampPalette(c("navy","blue","white","orange","red"))(100))
dev.off()


LN_TPM_ShanXi=read.table("D:/06_CRC/bulkRNAseq/ShanXiSample/LN_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
LN_Info_ShanXi=read.table("D:/06_CRC/bulkRNAseq/ShanXiSample/SubTypeIdentification/LN_SubType_Info.txt",header=T,row.names=1,sep="\t")
LN_Info_ShanXi$GroupByGene=factor(LN_Info_ShanXi$GroupByGene,levels=sort(unique(LN_Info_ShanXi$GroupByGene)))
LN_Info_ShanXi=LN_Info_ShanXi[order(LN_Info_ShanXi$GroupByGene),]
LN_TPM_ShanXi=LN_TPM_ShanXi[,rownames(LN_Info_ShanXi)]
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info_ShanXi$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info_ShanXi$GroupByGene))

targetTF=c("EGR1","CEBPD","SOX4","HES1","KLF10","HIF1A","ELF3","PRDM1","ETS2","XBP1","CEBPB","FOXF2","NFIC","FOXF1","PRNP","SMARCA1","NFIX","STAT1","SPIB","LTF","RELB","MEF2C","KLF4","ERG","SMAD1","MECOM","SOX18")
targetTFExpr=LN_TPM_ShanXi[targetTF,]
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info_ShanXi$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info_ShanXi$GroupByGene))
ann_colors = list(
                    GroupByGene = GroupByGene_Colors, 
                    LNGroup=c("Pos"="Orange","Neg"="LightBlue")
                )
pdf("D:/06_CRC/Graph/Part4/TFExprInShanXi.pdf",width=8,height=4)
pheatmap(targetTFExpr,annotation_col=LN_Info_ShanXi[,c("GroupByGene","LNGroup")],gaps_row = c(17,21),
    cluster_col=F,cluster_row=F,clustering_method="ward.D2",scale="row",annotation_colors = ann_colors,
    show_colnames=F,color = colorRampPalette(c("navy","blue","white","orange","red"))(100))
dev.off()


LN_TPM=read.table("D:/06_CRC/bulkRNAseq/HarBinSample/LNSample_Gene_TPM_Syumbol.txt",header=T,row.names=1,sep="\t")
LN_Info=read.table("D:/06_CRC/bulkRNAseq/HarBinSample/LN_GroupInfo.txt",header=T,row.names=1,sep="\t")
LN_Info$GroupByGene=factor(LN_Info$GroupByGene,levels=sort(unique(LN_Info$GroupByGene)))
LN_Info=LN_Info[order(LN_Info$GroupByGene),]
LN_TPM=LN_TPM[,rownames(LN_Info)]
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info$GroupByGene))

targetTF=c("EGR1","CEBPD","SOX4","HES1","KLF10","HIF1A","ELF3","PRDM1","ETS2","XBP1","CEBPB","FOXF2","NFIC","FOXF1","PRNP","SMARCA1","NFIX","STAT1","SPIB","LTF","RELB","MEF2C","KLF4","ERG","SMAD1","MECOM","SOX18")
targetTFExpr=LN_TPM[targetTF,]
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info$GroupByGene))
ann_colors = list(
                    GroupByGene = GroupByGene_Colors, 
                    LNGroup=c("Pos"="Orange","Neg"="LightBlue")
                )
pdf("D:/06_CRC/Graph/Part4/TFExprInHarBin.pdf",width=8,height=4)
pheatmap(targetTFExpr,annotation_col=LN_Info[,c("GroupByGene","LNGroup")],gaps_row = c(17,21),
	cluster_col=F,cluster_row=F,clustering_method="ward.D2",scale="row",annotation_colors = ann_colors,
	show_colnames=F,color = colorRampPalette(c("navy","blue","white","orange","red"))(100))
dev.off()







# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="SPIB"]

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
SPIBTableSubset <- regulonTargetsInfo[TF=="SPIB" & highConfAnnot==TRUE]
head(SPIBTableSubset)
write.table(SPIBTableSubset$gene,file="/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/SPIB.targetGene.txt",quote=F,row.names=F,col.names=F,sep="\t")

RELBTableSubset <- regulonTargetsInfo[TF=="RELB" & highConfAnnot==TRUE]
sort(table(RELBTableSubset$bestMotif),decreasing=TRUE)
write.table(RELBTableSubset$gene,file="/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/RELB.targetGene.txt",quote=F,row.names=F,col.names=F,sep="\t")

LTFTableSubset <- regulonTargetsInfo[TF=="LTF" & highConfAnnot==TRUE]
sort(table(LTFTableSubset$bestMotif),decreasing=TRUE)
write.table(LTFTableSubset$gene,file="/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/LTF.targetGene.txt",quote=F,row.names=F,col.names=F,sep="\t")













# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="SPIB"]
viewMotifs(tableSubset) 

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

# Cell-type specific regulators (RSS): 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)








##利用Seurat可视化AUC
dir.create('scenic_seurat')
#FeaturePlot
p1 = FeaturePlot(scRNAauc, features='CEBPB_extended_2290g', label=T, reduction = 'tsne')
p2 = FeaturePlot(scRNAbin, features='CEBPB_extended_2290g', label=T, reduction = 'tsne')
p3 = DimPlot(scRNA, reduction = 'tsne', group.by = "celltype_Monaco", label=T)
plotc = p1|p2|p3
ggsave('scenic_seurat/CEBPB_extended_2290g.png', plotc, width=14 ,height=4)






subSource=subset(ExAll,cells = colnames(regulonAUC)) #for ExAll.integrated
selected_f_10xGenomic <- rownames(subSource)[Matrix::rowSums(subSource) > ncol(subSource)*0.1] #expressed in at lest 10% cells
print(paste0(s,":Expressed Genes ",length(selected_f_10xGenomic)))
#7665
cellInfo=data.frame(clusters=subSource$seurat_clusters,Layers=subSource$Layer,Gender=subSource$Gender,Age=subSource$Age,Status=subSource$Statues)
clusterNumber=table(cellInfo$Layers)
cellInfo=cellInfo[cellInfo$Layers%in%names(clusterNumber[clusterNumber>10]),]
#https://www.jianshu.com/p/7ab2d6c8f764
print(paste0(s,":TFNumber ",dim(regulonAUC)[1]))
#507 10000
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
cellId=intersect(rownames(cellInfo),colnames(regulonAUC))
regulonAUC <- regulonAUC[,cellId]
cellInfo=cellInfo[cellId,]
all(rownames(cellInfo)==colnames(regulonAUC))
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Layers), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType <- regulonActivity_byCellType[,unique(cellInfo$Layers)]
dim(regulonActivity_byCellType)
regulonActivity_byCellType=regulonActivity_byCellType[Biobase::rowMax(regulonActivity_byCellType)>0.01,] #only focused on the TFs with auc socre greater than 0.01 in at least one sub cluster
regulonActivity_byCellType.df=data.frame(regulonActivity_byCellType)
regulonActivity_byCellType.df$TFName=stringr::str_sub(rownames(regulonActivity_byCellType.df), end=-4)
regulonActivity_byCellType.df=regulonActivity_byCellType.df[regulonActivity_byCellType.df$TF%in%selected_f_10xGenomic,]
regulonActivity_byCellType=data.frame(regulonActivity_byCellType,check.names=F)
regulonActivity_byCellType=regulonActivity_byCellType[rownames(regulonActivity_byCellType.df),]
print(paste0(s,":TFNumber_AfterFilter: ",dim(regulonActivity_byCellType)[1]))
dim(regulonActivity_byCellType)
pdf(paste0("Heatmap/",s,"_FullName.pdf",sep=""),height=20)
t=pheatmap(regulonActivity_byCellType,scale="row",clustering_method="ward.D2")
dev.off()














library(Seurat)
library(SeuratDisk)
library(SeuratData)


EndAndFib.loom <- as.loom(EndAndFib, filename = "/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC/EndAndFib.loom", verbose = FALSE)
EndAndFib.loom$close_all()


pyscenic grn \
${source}.downsampling.loom \
--num_workers 5 \
--output ${source}.downsampling.tsv \
--method grnboost2 \
/Projects/deng/Public/SCENIC/allTFs_hg38.txt

pyscenic ctx \
${source}.downsampling.tsv \
/Projects/deng/Public/SCENIC/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather /Projects/deng/Public/SCENIC/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname /Projects/deng/Aging/Ex/AllEx/SCENIC/databases/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname ${source}.downsampling.loom \
--transpose \
--mask_dropouts \
--mode "dask_multiprocessing" \
--output ${source}.reg.csv \
--num_workers 5

pyscenic aucell \
${source}.downsampling.loom \
${source}.reg.csv \
--output ${source}.downsampling.final.loom \
--num_workers 5




setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/scRNAseq/EndAndFib/SCENIC")
### Initialize settings
org <- "hgnc" # or hgnc, or dmel
dbDir <- "/lkn_lab/bixiaoman/Software/SCENIC/cisTarget_databases"
myDatasetTitle <- "SCENIC" # choose a name for your analysis
data(defaultDbNames)
#https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/
defaultDbNames$hgnc["500bp"] <- "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
defaultDbNames$hgnc["10kb"] <- "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
defaultDbNames
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
scenicOptions <- initializeScenic(org = org,
                                  dbDir = dbDir,
                                  dbs = defaultDbNames[["hgnc"]],
                                  datasetTitle = myDatasetTitle,
                                  nCores = 10)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$seed <- 123
