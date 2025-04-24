library(ggalluvial)
library(ggplot2)
#########  PCA analyis of both LN and Tissue  ##########################
GeneTPM_Final_Symbol=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/AllSample_Gene_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
AllSampleInfo=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SampleInformation_Final.txt",sep="\t",header=T,row.names=1)
IDMapping=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SampleAndMedicalID.txt",header=T,sep="\t")
rownames(IDMapping)=IDMapping$Sample.ID

TissueInfo=AllSampleInfo[AllSampleInfo$Tissue%in%c("Tumor"),]
TissueSampleList=intersect(rownames(TissueInfo),colnames(GeneTPM_Final_Symbol))
length(TissueSampleList)
#85

TissueSampleWithMedicalRecord=intersect(TissueSampleList,rownames(IDMapping))
length(TissueSampleWithMedicalRecord)
#85
IDMapping=IDMapping[TissueSampleWithMedicalRecord,]
TissueGeneExpr=GeneTPM_Final_Symbol[,TissueSampleWithMedicalRecord]
all(rownames(IDMapping)==colnames(TissueGeneExpr))
#TRUE
colnames(TissueGeneExpr)=IDMapping$MedicalRecordID
write.table(TissueGeneExpr,file="D:/06_CRC/bulkRNAseq/BeiJingSample/MSS/PrimaryTissueGeneExpr.txt",quote=F,sep="\t")


library(CMScaller)
Beijing_Tissue_TPM=read.table("/boot3/bixm/CRCTmp/PrimaryTissueGeneExpr.txt",header=T,row.names=1,sep="\t",check.names=F)
res=CMScaller(
       Beijing_Tissue_TPM,
       templates = CMScaller::templates.CMS,
       rowNames = "symbol",
       RNAseq = TRUE,
       nPerm = 1000,
       seed = NULL,
       FDR = 0.05,
       doPlot = TRUE,
       verbose = TRUE
     )
#CMS1 CMS2 CMS3 CMS4 <NA> 
#   6   15   14   20    9 
#9/64 samples set to NA
write.table(res,file="/boot3/bixm/CRCTmp/CMSPredication4PrimaryTissue.txt",sep="\t",quote=F)

BeiJingInfo=read.table("/boot3/bixm/CRCTmp/LymphNodeInfomation-Beijing_173.txt",header=T,sep="\t")
BeiJingInfo$CombineGroup=ifelse(BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2"),"NLN_C1/C2",ifelse(BeiJingInfo$GroupByGene%in%c("NLN_C3","NLN_C4"),"NLN_C3/C4","PLN"))
TargetInfo=c("MedicalRecordID","CombineGroup","GroupByGene","LNSize","LNLocation","TumorLocation","Age","Sex","BMI","DifferentiationGrade","PathologicalType","TumorSize","TNM.stage","MSS.MSI","BRAF","Vascular.Thrombosis","PNI","EMVI","Tumor.Budding","OS","RFS")
BeiJingInfo_SubType=BeiJingInfo[BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2","NLN_C3","NLN_C4"),TargetInfo]

CMSResult=read.table("/boot3/bixm/CRCTmp/CMSPredication4PrimaryTissue.txt",header=T,sep="\t")
CMSResult$MedicalRecordID=rownames(CMSResult)
BeiJingInfo_SubType=merge(BeiJingInfo_SubType,CMSResult,by="MedicalRecordID",all.x=TRUE)

GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246")
names(GroupByGene_Colors)=sort(unique(BeiJingInfo_SubType$GroupByGene))
table(BeiJingInfo_SubType$prediction,BeiJingInfo_SubType$GroupByGene)
BeiJingInfo_SubType[is.na(BeiJingInfo_SubType)]<-"Not defined"
RelationData=data.frame(table(BeiJingInfo_SubType$prediction,BeiJingInfo_SubType$GroupByGene))
colnames(RelationData)=c("CMSSubtype","SubType","Number")

g=ggplot(data = RelationData,aes(axis1 = SubType, axis2 = CMSSubtype, y = Number)) +
  geom_alluvium(aes(fill = SubType)) +
  geom_stratum()+
  scale_fill_manual(values=GroupByGene_Colors) +
  geom_text(stat = "stratum",aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("CMSSubtype", "SubType"),expand = c(0.15, 0.05)) +
  theme_void()
pdf("/boot3/bixm/CRCTmp/SubTypeAndCMS_PrimaryTissue.pdf",width=6,height=5)
print(g)
dev.off()

GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246")
names(GroupByGene_Colors)=sort(unique(BeiJingInfo_SubType$GroupByGene))
CMSResult=data.frame(table(BeiJingInfo_SubType$CombineGroup,BeiJingInfo_SubType$prediction))
colnames(CMSResult)=c("SubTypes","CMSResult","Number")
g=ggplot(CMSResult,aes(x=SubTypes, y=Number, fill=CMSResult)) +
      geom_bar(position = 'fill',stat = "identity") +
      #scale_fill_manual(values=GroupByGene_Colors) +
      scale_fill_brewer(palette = "Dark2")+
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("/boot3/bixm/CRCTmp/SubTypeAndCMS_PrimaryTissue.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
CMSResult.mt=as.data.frame.matrix(table(BeiJingInfo_SubType$CombineGroup,BeiJingInfo_SubType$prediction))
fisher.test(CMSResult.mt) #p-value = 0.1494


GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246")
names(GroupByGene_Colors)=sort(unique(BeiJingInfo_SubType$GroupByGene))

BeiJingInfo_SubType=BeiJingInfo_SubType[BeiJingInfo_SubType$prediction%in%c("CMS1","CMS2","CMS3","CMS4"),]
CMSResult=data.frame(table(BeiJingInfo_SubType$CombineGroup,BeiJingInfo_SubType$prediction))
colnames(CMSResult)=c("SubTypes","CMSResult","Number")
g=ggplot(CMSResult,aes(x=SubTypes, y=Number, fill=CMSResult)) +
      geom_bar(position = 'fill',stat = "identity") +
      #scale_fill_manual(values=GroupByGene_Colors) +
      scale_fill_brewer(palette = "Dark2")+
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("/boot3/bixm/CRCTmp/SubTypeAndCMS_PrimaryTissue.barplot.noNA.pdf",width=3,height=3)
print(g)
dev.off()
CMSResult.mt=as.data.frame.matrix(table(BeiJingInfo_SubType$CombineGroup,BeiJingInfo_SubType$prediction))
fisher.test(CMSResult.mt) #p-value = 0.1494






#intrinsic CMS (iCMS) subtypes?
library(CMScaller)
Beijing_Tissue_TPM=read.table("/boot3/bixm/CRCTmp/PrimaryTissueGeneExpr.txt",header=T,row.names=1,sep="\t",check.names=F)
iCMSMarker=read.table("/boot3/bixm/CRCTmp/iCMSMarkerForCMSCaller.txt",header=T,sep="\t")
head(iCMSMarker)
res=CMScaller(
       Beijing_Tissue_TPM,
       templates = iCMSMarker,
       rowNames = "symbol",
       RNAseq = TRUE,
       nPerm = 1000,
       seed = NULL,
       FDR = 0.05,
       doPlot = TRUE,
       verbose = TRUE
     )

#iCMS2 iCMS3  <NA> 
#   28    27     9
write.table(res,file="/boot3/bixm/CRCTmp/iCMSPredication4PrimaryTissue.txt",sep="\t",quote=F)


BeiJingInfo=read.table("/boot3/bixm/CRCTmp/LymphNodeInfomation-Beijing_173.txt",header=T,sep="\t")
BeiJingInfo$CombineGroup=ifelse(BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2"),"NLN_C1/C2",ifelse(BeiJingInfo$GroupByGene%in%c("NLN_C3","NLN_C4"),"NLN_C3/C4","PLN"))
TargetInfo=c("MedicalRecordID","CombineGroup","GroupByGene","LNSize","LNLocation","TumorLocation","Age","Sex","BMI","DifferentiationGrade","PathologicalType","TumorSize","TNM.stage","MSS.MSI","BRAF","Vascular.Thrombosis","PNI","EMVI","Tumor.Budding","OS","RFS")
BeiJingInfo_SubType=BeiJingInfo[BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2","NLN_C3","NLN_C4"),TargetInfo]

CMSResult=read.table("/boot3/bixm/CRCTmp/iCMSPredication4PrimaryTissue.txt",header=T,sep="\t")
CMSResult$MedicalRecordID=rownames(CMSResult)
BeiJingInfo_SubType=merge(BeiJingInfo_SubType,CMSResult,by="MedicalRecordID",all.x=TRUE)

GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246")
names(GroupByGene_Colors)=sort(unique(BeiJingInfo_SubType$GroupByGene))
table(BeiJingInfo_SubType$prediction,BeiJingInfo_SubType$GroupByGene)
BeiJingInfo_SubType[is.na(BeiJingInfo_SubType)]<-"Not defined"
RelationData=data.frame(table(BeiJingInfo_SubType$prediction,BeiJingInfo_SubType$GroupByGene))
colnames(RelationData)=c("CMSSubtype","SubType","Number")

g=ggplot(data = RelationData,aes(axis1 = SubType, axis2 = CMSSubtype, y = Number)) +
  geom_alluvium(aes(fill = SubType)) +
  geom_stratum()+
  scale_fill_manual(values=GroupByGene_Colors) +
  geom_text(stat = "stratum",aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("CMSSubtype", "SubType"),expand = c(0.15, 0.05)) +
  theme_void()
pdf("/boot3/bixm/CRCTmp/SubTypeAndiCMS_PrimaryTissue.pdf",width=6,height=5)
print(g)
dev.off()

GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246")
names(GroupByGene_Colors)=sort(unique(BeiJingInfo_SubType$GroupByGene))
CMSResult=data.frame(table(BeiJingInfo_SubType$CombineGroup,BeiJingInfo_SubType$prediction))
colnames(CMSResult)=c("SubTypes","CMSResult","Number")
g=ggplot(CMSResult,aes(x=SubTypes, y=Number, fill=CMSResult)) +
      geom_bar(position = 'fill',stat = "identity") +
      #scale_fill_manual(values=GroupByGene_Colors) +
      scale_fill_brewer(palette = "Dark2")+
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("/boot3/bixm/CRCTmp/SubTypeAndiCMS_PrimaryTissue.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
CMSResult.mt=as.data.frame.matrix(table(BeiJingInfo_SubType$CombineGroup,BeiJingInfo_SubType$prediction))
fisher.test(CMSResult.mt) #p-value = 0.1494

BeiJingInfo_SubType=BeiJingInfo_SubType[BeiJingInfo_SubType$prediction%in%c("iCMS2","iCMS3"),]
CMSResult=data.frame(table(BeiJingInfo_SubType$CombineGroup,BeiJingInfo_SubType$prediction))
colnames(CMSResult)=c("SubTypes","CMSResult","Number")
g=ggplot(CMSResult,aes(x=SubTypes, y=Number, fill=CMSResult)) +
      geom_bar(position = 'fill',stat = "identity") +
      #scale_fill_manual(values=GroupByGene_Colors) +
      scale_fill_brewer(palette = "Dark2")+
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("/boot3/bixm/CRCTmp/SubTypeAndiCMS_PrimaryTissue.barplot.noNA.pdf",width=3,height=3)
print(g)
dev.off()
CMSResult.mt=as.data.frame.matrix(table(BeiJingInfo_SubType$CombineGroup,BeiJingInfo_SubType$prediction))
fisher.test(CMSResult.mt) #p-value = 0.1494
















#Does the distribution of non-metastatic lymph nodes correlate with Consensus Molecular Subtypes (CMS) or intrinsic CMS (iCMS) subtypes?
#Consensus Molecular Subtypes (CMS)
library(CMScaller)
Beijing_LN_TPM=read.table("/boot3/bixm/CRCTmp/LN_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
res=CMScaller(
       Beijing_LN_TPM,
       templates = CMScaller::templates.CMS,
       rowNames = "symbol",
       RNAseq = TRUE,
       nPerm = 1000,
       seed = NULL,
       FDR = 0.05,
       doPlot = TRUE,
       verbose = TRUE
     )
#CMS1 CMS2 CMS3 CMS4 <NA> 
#  6   14    8   58   87 
#87/173 samples set to NA
write.table(res,file="/boot3/bixm/CRCTmp/CMSPredication.txt",sep="\t",quote=F)

res=read.table("/boot3/bixm/CRCTmp/CMSPredication.txt",header=T,sep="\t")
LN_Info=read.table("/boot3/bixm/CRCTmp/LN_SubType_Info.txt",header=T,row.names=1,sep="\t")
LN_CMS=res[rownames(LN_Info),]
all(rownames(LN_CMS)==rownames(LN_Info))
#TRUE
LN_CMS_Info=data.frame(cbind(LNId=rownames(LN_CMS),CMSSubtype=LN_CMS$prediction,CMSpvalue=LN_CMS$p.value,SubType=LN_Info$GroupByGene))

GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246","#E2BF8C","#147286")
names(GroupByGene_Colors)=sort(unique(LN_CMS_Info$GroupByGene))
table(LN_CMS_Info$CMSSubtype,LN_CMS_Info$SubType)
LN_CMS_Info[is.na(LN_CMS_Info)]<-"Not defined"
RelationData=data.frame(table(LN_CMS_Info$CMSSubtype,LN_CMS_Info$SubType))
colnames(RelationData)=c("CMSSubtype","SubType","Number")
RelationData=RelationData[RelationData$SubType%in%c("NLN_C1","NLN_C2","NLN_C3","NLN_C4"),]
g=ggplot(data = RelationData,aes(axis1 = SubType, axis2 = CMSSubtype, y = Number)) +
  geom_alluvium(aes(fill = SubType)) +
  geom_stratum()+
  scale_fill_manual(values=GroupByGene_Colors) +
  geom_text(stat = "stratum",aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("CMSSubtype", "SubType"),expand = c(0.15, 0.05)) +
  theme_void()
pdf("/boot3/bixm/CRCTmp/SubTypeAndCMS.pdf",width=6,height=5)
print(g)
dev.off()

#intrinsic CMS (iCMS) subtypes?
library(CMScaller)
iCMSMarker=read.table("/boot3/bixm/CRCTmp/iCMSMarkerForCMSCaller.txt",header=T,sep="\t")
head(iCMSMarker)
Beijing_LN_TPM=read.table("/boot3/bixm/CRCTmp/LN_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
res=CMScaller(
       Beijing_LN_TPM,
       templates = iCMSMarker,
       rowNames = "symbol",
       RNAseq = TRUE,
       nPerm = 1000,
       seed = NULL,
       FDR = 0.05,
       doPlot = TRUE,
       verbose = TRUE
     )

#iCMS2 iCMS3  <NA> 
#  53    55    65
#87/173 samples set to NA
write.table(res,file="/boot3/bixm/CRCTmp/iCMSPredication.txt",sep="\t",quote=F)

res=read.table("/boot3/bixm/CRCTmp/iCMSPredication.txt",header=T,sep="\t")
LN_Info=read.table("/boot3/bixm/CRCTmp/LN_SubType_Info.txt",header=T,row.names=1,sep="\t")
LN_CMS=res[rownames(LN_Info),]
all(rownames(LN_CMS)==rownames(LN_Info))
#TRUE
LN_CMS_Info=data.frame(cbind(LNId=rownames(LN_CMS),CMSSubtype=LN_CMS$prediction,CMSpvalue=LN_CMS$p.value,SubType=LN_Info$GroupByGene))

GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246","#E2BF8C","#147286")
names(GroupByGene_Colors)=sort(unique(LN_CMS_Info$GroupByGene))
table(LN_CMS_Info$CMSSubtype,LN_CMS_Info$SubType)
LN_CMS_Info[is.na(LN_CMS_Info)]<-"Not defined"
RelationData=data.frame(table(LN_CMS_Info$CMSSubtype,LN_CMS_Info$SubType))
colnames(RelationData)=c("CMSSubtype","SubType","Number")
RelationData=RelationData[RelationData$SubType%in%c("NLN_C1","NLN_C2","NLN_C3","NLN_C4"),]
g=ggplot(data = RelationData,aes(axis1 = SubType, axis2 = CMSSubtype, y = Number)) +
  geom_alluvium(aes(fill = SubType)) +
  geom_stratum()+
  scale_fill_manual(values=GroupByGene_Colors) +
  geom_text(stat = "stratum",aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("CMSSubtype", "SubType"),expand = c(0.15, 0.05)) +
  theme_void()
pdf("/boot3/bixm/CRCTmp/SubTypeAndiCMS.pdf",width=6,height=5)
print(g)
dev.off()

