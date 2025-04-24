library(Seurat)
library(monocle)
library(RColorBrewer)

setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample")
Beijing_LN_Count=read.table("LN_Count_Final_Symbol.txt",header=T,row.names=1,sep="\t")
Beijing_LN_Info=read.table("LN_SubType_Info.txt",header=T,row.names=1,sep="\t")
Beijing_LN_Count=Beijing_LN_Count[,rownames(Beijing_LN_Info)]
all(rownames(Beijing_LN_Info)==colnames(Beijing_LN_Count))

setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/ShanXiSample")
Shanxi_LN_Count=read.table("LNSample_Gene_Count_Symbol.txt",header=T,row.names=1,sep="\t")
Shanxi_LN_Info=read.table("LN_GroupInfo.txt",header=T,row.names=1,sep="\t")
Shanxi_LN_Count=Shanxi_LN_Count[,rownames(Shanxi_LN_Info)]
all(rownames(Shanxi_LN_Info)==colnames(Shanxi_LN_Count))

setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/HarBinSample")
Harbin_LN_Count=read.table("LNSample_Gene_Count_Symbol.txt",header=T,row.names=1,sep="\t")
Harbin_LN_Info=read.table("LN_GroupInfo.txt",header=T,row.names=1,sep="\t")
Harbin_LN_Count=Harbin_LN_Count[,rownames(Harbin_LN_Info)]
all(rownames(Harbin_LN_Info)==colnames(Harbin_LN_Count))

Beijing_LN_Info$Dataset="Beijing"
Shanxi_LN_Info$Dataset="Shanxi"
Harbin_LN_Info$Dataset="Harbin"
Beijing_LN_Info$LNStatus=Beijing_LN_Info$LNGroup
AllLNInfo=rbind(Beijing_LN_Info[,c("GroupByGene","LNStatus","Dataset")],Shanxi_LN_Info[,c("GroupByGene","LNStatus","Dataset")],Harbin_LN_Info[,c("GroupByGene","LNStatus","Dataset")])

GeneList=Reduce(intersect,list(rownames(Beijing_LN_Count),rownames(Shanxi_LN_Count),rownames(Harbin_LN_Count)))
length(GeneList)
#56779

dim(Beijing_LN_Count)
dim(Shanxi_LN_Count)
dim(Harbin_LN_Count)
AllLNCount=cbind(Beijing_LN_Count[GeneList,],Shanxi_LN_Count[GeneList,],Harbin_LN_Count[GeneList,])
all(rownames(AllLNInfo)==colnames(AllLNCount))

logTPM=log2(AllLNCount+1)
logTPM=logTPM[rowMeans(logTPM)>0,]
logTPM.pca <- prcomp(t(logTPM))
head(logTPM.pca$x)
all(rownames(logTPM.pca$x)==rownames(AllLNInfo))
df_pca=data.frame(PC1=logTPM.pca$x[,1],PC2=logTPM.pca$x[,2],GroupByGene=AllLNInfo$GroupByGene,LNStatus=AllLNInfo$LNStatus,Dataset=AllLNInfo$Dataset)
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(AllLNInfo$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(AllLNInfo$GroupByGene))
g=ggplot(df_pca,aes(x=PC1,y=PC2,color=Dataset,shape=LNStatus))+ 
geom_point(size=3)+
theme_bw()+
scale_shape_manual(values=c(19, 8))
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/Monocle4All/LN_PCA.pdf",width=4,height=4)
print(g)
dev.off()

mod=model.matrix(~factor(AllLNInfo$GroupByGene))
AllLNCount_BatchRemoved=limma::removeBatchEffect(AllLNCount, AllLNInfo$Dataset,design=mod)
AllLNCount_BatchRemoved.integer=round(AllLNCount_BatchRemoved)
AllLNCount_BatchRemoved.integer[AllLNCount_BatchRemoved.integer<0]<-0
logTPM=log2(AllLNCount_BatchRemoved.integer+1)
logTPM=logTPM[rowMeans(logTPM)>0,]
logTPM.pca <- prcomp(t(logTPM))
head(logTPM.pca$x)
all(rownames(logTPM.pca$x)==rownames(AllLNInfo))
df_pca=data.frame(PC1=logTPM.pca$x[,1],PC2=logTPM.pca$x[,2],GroupByGene=AllLNInfo$GroupByGene,LNStatus=AllLNInfo$LNStatus,Dataset=AllLNInfo$Dataset)
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(AllLNInfo$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(AllLNInfo$GroupByGene))
g=ggplot(df_pca,aes(x=PC1,y=PC2,color=Dataset,shape=LNStatus))+ 
geom_point(size=3)+
theme_bw()+
scale_shape_manual(values=c(19, 8))
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/Monocle4All/LN_PCA_BatchRemoved.pdf",width=5,height=4)
print(g)
dev.off()


AllLNInfo=rbind(Shanxi_LN_Info[,c("GroupByGene","LNStatus","Dataset")],Harbin_LN_Info[,c("GroupByGene","LNStatus","Dataset")])
GeneList=Reduce(intersect,list(rownames(Shanxi_LN_Count),rownames(Harbin_LN_Info)))
length(GeneList)
#56779
dim(Beijing_LN_Count)
dim(Shanxi_LN_Count)
dim(Harbin_LN_Count)
AllLNCount=cbind(Shanxi_LN_Count[GeneList,],Harbin_LN_Count[GeneList,])
all(rownames(AllLNInfo)==colnames(AllLNCount))


all(rownames(Harbin_LN_Info)==colnames(Harbin_LN_Count))
p_data <- Harbin_LN_Info
f_data <- data.frame(gene_short_name=rownames(Harbin_LN_Count),row.names=rownames(Harbin_LN_Count))
pd <- new("AnnotatedDataFrame",data=p_data)
fd <- new("AnnotatedDataFrame",data=f_data)
cds <- newCellDataSet(as.matrix(Harbin_LN_Count),phenoData=pd,featureData=fd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
length(expressed_genes)
#45431
diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~GroupByGene")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds,reverse =TRUE)


GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Dark2"))(length(unique(AllLNInfo$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(AllLNInfo$GroupByGene))
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/Monocle4All/Harbin_LN_SubType.pdf",width=4,height=4.3)
plot_cell_trajectory(cds, color_by = "GroupByGene")&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/Monocle4All/Harbin_LN_Pseudotime.pdf",width=4,height=4)
plot_cell_trajectory(cds, color_by = "Pseudotime")&scale_color_viridis()&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/Monocle4All/Harbin_LN_Gene_MKI67.pdf",width=4,height=4)
plot_cell_trajectory(cds, markers="MKI67")&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/Monocle4All/Harbin_LN_Dataset.pdf",width=4,height=4.3)
plot_cell_trajectory(cds, color_by = "Dataset")&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()

saveRDS(cds,file="Harbin_CohortCombine.rds")







all(rownames(Shanxi_LN_Info)==colnames(Shanxi_LN_Count))
p_data <- Shanxi_LN_Info
f_data <- data.frame(gene_short_name=rownames(Shanxi_LN_Count),row.names=rownames(Shanxi_LN_Count))
pd <- new("AnnotatedDataFrame",data=p_data)
fd <- new("AnnotatedDataFrame",data=f_data)
cds <- newCellDataSet(as.matrix(Shanxi_LN_Count),phenoData=pd,featureData=fd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
length(expressed_genes)
#45431
diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~GroupByGene")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds,reverse =TRUE)


GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Dark2"))(length(unique(AllLNInfo$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(AllLNInfo$GroupByGene))
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/Monocle4All/Shanxi_LN_SubType.pdf",width=4,height=4.3)
plot_cell_trajectory(cds, color_by = "GroupByGene")&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/Monocle4All/Shanxi_LN_Pseudotime.pdf",width=4,height=4)
plot_cell_trajectory(cds, color_by = "Pseudotime")&scale_color_viridis()&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/Monocle4All/Shanxi_LN_Gene_MKI67.pdf",width=4,height=4)
plot_cell_trajectory(cds, markers="MKI67")&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()
pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/Monocle4All/Shanxi_LN_Dataset.pdf",width=4,height=4.3)
plot_cell_trajectory(cds, color_by = "Dataset")&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()

saveRDS(cds,file="Shanxi_CohortCombine.rds")


