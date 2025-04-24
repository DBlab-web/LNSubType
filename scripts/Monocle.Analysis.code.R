salloc -N 1 -n 32 -p 128c1t --comment=lkn_lab
source /opt/app/anaconda3/bin/activate
conda activate Seurat-4.4.0

library(Seurat)
library(monocle)
library(RColorBrewer)

setwd("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample")
LN_Count=read.table("LN_Count_Final_Symbol.txt",header=T,row.names=1,sep="\t")
LN_Info=read.table("LN_SubType_Info.txt",header=T,row.names=1,sep="\t")
LN_Count=LN_Count[,rownames(LN_Info)]
all(rownames(LN_Info)==colnames(LN_Count))

p_data <- LN_Info
f_data <- data.frame(gene_short_name=rownames(LN_Count),row.names=rownames(LN_Count))

pd <- new("AnnotatedDataFrame",data=p_data)
fd <- new("AnnotatedDataFrame",data=f_data)

cds <- newCellDataSet(as.matrix(LN_Count),phenoData=pd,featureData=fd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))

cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))

diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~GroupByGene")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds,reverse =TRUE)

GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info$GroupByGene))

pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample/Monocle/LN_SubType.pdf",width=4,height=4.3)
plot_cell_trajectory(cds, color_by = "GroupByGene")&scale_color_manual(values=GroupByGene_Colors)&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()

pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample/Monocle/LN_Pseudotime.pdf",width=4,height=4)
plot_cell_trajectory(cds, color_by = "Pseudotime")&scale_color_viridis()&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()

pdf("/groups/g900008/home/bixiaoman/Projects/LNCRC/bulkRNAseq/BeiJingSample/Monocle/LN_Gene_MKI67.pdf",width=4,height=4)
plot_cell_trajectory(cds, markers="MKI67")&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.0, linetype="solid"))
dev.off()

saveRDS(cds,file="Monocle/LN_Pseudotime.rds")



