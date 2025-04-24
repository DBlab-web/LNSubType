
GeneInfo=read.csv("D:/06_CRC/bulkRNAseq/BeiJingSample/GeneInformation.txt",sep="\t",header=T,row.names=1)
protein_coding_genes=GeneInfo[GeneInfo$gene_biotype=="protein_coding","gene_name"]


PosvsNeg_DEG=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_PosvsNeg.txt",header=T,row.name=1,sep="\t")
NLNC4vsC1_DEG=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/AllGene_NLN_C4vsNLN_C1.txt",header=T,row.name=1,sep="\t")

PosvsNeg_DEG$threshold=ifelse(PosvsNeg_DEG$padj>0.01|abs(PosvsNeg_DEG$log2FoldChange)<log2(2),"NoSig",ifelse(PosvsNeg_DEG$log2FoldChange>log2(2),"UpInPos","DownInPos"))
NLNC4vsC1_DEG$threshold=ifelse(NLNC4vsC1_DEG$padj>0.01|abs(NLNC4vsC1_DEG$log2FoldChange)<log2(2),"NoSig",ifelse(NLNC4vsC1_DEG$log2FoldChange>log2(2),"UpInNLNC4","DownInNLNC4"))

PosvsNeg_DEG=PosvsNeg_DEG[intersect(rownames(PosvsNeg_DEG),protein_coding_genes),]
NLNC4vsC1_DEG=NLNC4vsC1_DEG[intersect(rownames(NLNC4vsC1_DEG),protein_coding_genes),]


sigPos=PosvsNeg_DEG[!PosvsNeg_DEG$threshold%in%"NoSig",]
sigC4=NLNC4vsC1_DEG[!NLNC4vsC1_DEG$threshold%in%"NoSig",]
library(ggvenn)
sigPosList=split(rownames(sigPos),sigPos$threshold)
sigC4List=split(rownames(sigC4),sigC4$threshold)
DEGList=c(sigPosList,sigC4List)

pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/DEGOverlap.pdf",width=5,height=5)
ggvenn(DEGList,c("UpInPos","UpInNLNC4","DownInNLNC4","DownInPos"),fill_color = c("orange", "red", "blue", "lightblue"))
dev.off()


UpInNLNC4Spe=setdiff(sigC4List$UpInNLNC4,rownames(sigPos))
UpInNLNC4Spe=UpInNLNC4Spe[1:2999]
DownInNLNC4Spe=setdiff(sigC4List$DownInNLNC4,rownames(sigPos))
UpInPosSep=setdiff(sigPosList$UpInPos,rownames(sigC4))
DownInPosSep=setdiff(sigPosList$DownInPos,rownames(sigC4))
SharedUp=intersect(sigC4List$UpInNLNC4,sigPosList$UpInPos)
SharedDown=intersect(sigC4List$DownInNLNC4,sigPosList$DownInPos)

DEGList=rbind(toString(UpInNLNC4Spe),toString(UpInPosSep),toString(SharedUp),
              toString(DownInNLNC4Spe),toString(DownInPosSep),toString(SharedDown))
rownames(DEGList)=c("UpInNLNC4","UpInPos","UpInBoth","DownInNLNC4","DownInPos","DownInBoth")
write.table(DEGList[c(1,2,3),],file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/metascape/NLNC4AndPos.Up.combined.txt",col.names=F,quote=F,sep="\t")
write.table(DEGList[c(4,5,6),],file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/metascape/NLNC4AndPos.Down.combined.txt",col.names=F,quote=F,sep="\t")



write.table(SharedUp,file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/metascape/NLNC4AndPos.SahredUp.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(SharedDown,file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/metascape/NLNC4AndPos.SahredDown.txt",row.names=F,col.names=F,quote=F,sep="\t")

write.table(UpInNLNC4Spe,file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/metascape/NLNC4AndPos.UpInNLNC4SpeTop3000.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(DownInNLNC4Spe,file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/metascape/NLNC4AndPos.DownInNLNC4Sep.txt",row.names=F,col.names=F,quote=F,sep="\t")

write.table(UpInPosSep,file="D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/metascape/NLNC4AndPos.UpInPos4Spe.txt",row.names=F,col.names=F,quote=F,sep="\t")



data=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/metascape/NLNC4AndPos.UpInPos4Spe/NLNC4AndPos.UpInPos4Spe.Go.txt",header=T,sep="\t")
data=data[data$GroupID%in%grep("Summary",data$GroupID,value=TRUE),]
colorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(pbmc.merge$orig.ident)))
data$logP=-data$LogP
data$Description=factor(data$Description,levels=rev(data$Description))
g=ggplot(data[c(1:5),], aes(Description, logP, fill=logP)) +
  geom_bar(stat="identity") +
  coord_flip()+
  theme_bw()+
  scale_fill_viridis()+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.y = element_blank())

pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/DEG/metascape/NLNC4AndPos.UpInPos4Spe.Go.pdf",width=6.5,height=2)
print(g)
dev.off()














library(ggrepel)
hs_data=data.frame(NLNC4vsC1_DEG)
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


