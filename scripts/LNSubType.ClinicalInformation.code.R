library(ggplot2)
library(RColorBrewer)

BeiJingInfo=read.table("D:/06_CRC/bulkRNAseq/LNInfo/LymphNodeInfomation-Beijing_173.txt",header=T,sep="\t")
BeiJingInfo$CombineGroup=ifelse(BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2"),"NLN_C1/C2",ifelse(BeiJingInfo$GroupByGene%in%c("NLN_C3","NLN_C4"),"NLN_C3/C4","PLN"))

##The authors should comment on the distribution of non-metastatic lymph nodes both across patients and within individual patients. 
##This would provide insight into whether their representation is uniform or biased. 
head(BeiJingInfo)
NLNSample=BeiJingInfo[BeiJingInfo$LNMStatus=="NLN",]
NLNSampleInEachPerson=data.frame(table(NLNSample$MedicalRecordID))
colnames(NLNSampleInEachPerson)=c("SampleID","NLNNumber")
NLNNumberCount=data.frame(table(NLNSampleInEachPerson$NLNNumber))
colnames(NLNNumberCount)=c("non-metastatic_number","PersonCount")
# NLNNumberCount
#  non-metastatic_number PersonCount
#1                     1          27
#2                     2          28
#3                     3           8
#4                     4           2
#5                     5           2
#6                     7           2

g=ggplot(NLNNumberCount,aes(x=`non-metastatic_number`, y=PersonCount, fill=`non-metastatic_number`)) +
      geom_bar(stat = "identity") +
      scale_fill_brewer(palette = "Set2") +
      #scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(color="black"))
pdf("D:/06_CRC/Graph/Part6/BeiJing_non-metastatic_SampleDis_barplot.pdf",width=5,height=2.5)
print(g)
dev.off()

GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246","#E2BF8C","#147286")
names(GroupByGene_Colors)=sort(unique(BeiJingInfo$GroupByGene))

N1=NLNSampleInEachPerson[NLNSampleInEachPerson$NLNNumber==1,]
N1Detail=NLNSample[NLNSample$MedicalRecordID%in%N1$SampleID,]
N1SubType=data.frame(table(N1Detail$GroupByGene))
colnames(N1SubType)=c("SubTypes","Number")
pie <- ggplot(N1SubType, aes(x = "", y = Number, fill = SubTypes)) +
 geom_col(width = 1) +
 scale_fill_manual(values = c(GroupByGene_Colors)) +
 coord_polar(theta = "y")+
 theme_bw()+
 theme(panel.spacing = unit(0.05, "cm", data = NULL))+
 theme(axis.title.x=element_blank(),axis.title.y=element_blank())
pdf("D:/06_CRC/Graph/Part6/BeiJing_non-metastatic_SampleDis_N1_PiePlot.pdf",width=4,height=2.5)
print(pie)
dev.off()

> N1SubType
  SubTypes Number
1   NLN_C1      6
2   NLN_C2      4
3   NLN_C3      7
4   NLN_C4      6
5   PLN_C1      3
6   PLN_C2      1

#  SubTypes Number
1   NLN_C1     22
2   NLN_C2     13
3   NLN_C3     12
4   NLN_C4      6
5   PLN_C1      2
6   PLN_C2      1

> N3SubType
  SubTypes Number
1   NLN_C1      7
2   NLN_C2      8
3   NLN_C3      4
4   NLN_C4      5

> N4SubType
  SubTypes Number
1   NLN_C1      3
2   NLN_C2      1
3   NLN_C3      1
4   NLN_C4      3

> N5SubType
  SubTypes Number
1   NLN_C1      5
2   NLN_C2      3
3   NLN_C4      2

> N7SubType
  SubTypes Number
1   NLN_C1      5
2   NLN_C2      3
3   NLN_C3      5
4   NLN_C4      1



##How do non-metastatic lymph nodes relate to MSI/MSS status
TargetInfo=c("MedicalRecordID","CombineGroup","GroupByGene","LNSize","LNLocation","TumorLocation","Age","Sex","BMI","DifferentiationGrade","PathologicalType","HistologicalType","TumorSize","TNM.stage","MSS.MSI","BRAF","Vascular.Thrombosis","PNI","EMVI","Tumor.Budding","OS","RFS")
BeiJingInfo_NLN=BeiJingInfo[BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2","NLN_C3","NLN_C4"),TargetInfo]
MSIType=data.frame(table(BeiJingInfo_NLN$GroupByGene,BeiJingInfo_NLN$MSS.MSI))

colnames(MSIType)=c("SubTypes","MSIType","Number")
MSIType
#  SubTypes MSIType Number
#1   NLN_C1    dMMR      2
#2   NLN_C2    dMMR      2
#3   NLN_C3    dMMR      3
#4   NLN_C4    dMMR      1
#5   NLN_C1    pMMR     48
#6   NLN_C2    pMMR     33
#7   NLN_C3    pMMR     28
#8   NLN_C4    pMMR     31
g=ggplot(MSIType,aes(x=MSIType, y=Number, fill=SubTypes)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_Non_Metastatics_MSIType.barplot.pdf",width=3,height=4)
print(g)
dev.off()
MSIType.mt=as.data.frame.matrix(table(BeiJingInfo_NLN$GroupByGene,BeiJingInfo_NLN$MSS.MSI))
fisher.test(MSIType.mt) #p-value = 0.6504

#and other molecular or histopathological features of the primary tumor?
HistologicalType=data.frame(table(BeiJingInfo_NLN$GroupByGene,BeiJingInfo_NLN$HistologicalType))
colnames(HistologicalType)=c("SubTypes","HistologicalType","Number")
HistologicalType
#  SubTypes HistologicalType Number
#1   NLN_C1    Polypoid type     30
#2   NLN_C2    Polypoid type     22
#3   NLN_C3    Polypoid type     11
#4   NLN_C4    Polypoid type     15
#5   NLN_C1  Ulcerative type     20
#6   NLN_C2  Ulcerative type     13
#7   NLN_C3  Ulcerative type     20
#8   NLN_C4  Ulcerative type     17
g=ggplot(HistologicalType,aes(x=HistologicalType, y=Number, fill=SubTypes)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_Non_Metastatics_HistologicalType.barplot.pdf",width=3,height=4.5)
print(g)
dev.off()
HistologicalType.mt=as.data.frame.matrix(table(BeiJingInfo_NLN$GroupByGene,BeiJingInfo_NLN$HistologicalType))
fisher.test(HistologicalType.mt) #p-value = 0.08623



###### Figure 6 heatmap #################
SubtypeInfo=as.data.frame.matrix(table(BeiJingInfo$MedicalRecordID,BeiJingInfo$GroupByGene))
TargetInfo=c("MedicalRecordID","Age","Sex","BMI","OS","RFS","TumorLocation","PathologicalType","HistologicalType","DifferentiationGrade","TumorSize","TNM.stage","MSS.MSI","BRAF","Vascular.Thrombosis","PNI","EMVI","Tumor.Budding")
BeijingSampleInfo=unique(BeiJingInfo[,rev(TargetInfo)])
dim(BeijingSampleInfo)
#78 18
length(unique(BeijingSampleInfo$MedicalRecordID))
t=table(BeijingSampleInfo$MedicalRecordID)
#check the sampel with multiple MedicalRecordID
BeijingSampleInfo[BeijingSampleInfo$MedicalRecordID%in%names(t[t==2]),]
#none
rownames(BeijingSampleInfo)=BeijingSampleInfo$MedicalRecordID
BeijingSampleInfo$MedicalRecordID=NULL
# Specify colors
ann_colors = list(
    TumorLocation = c("Left colon" = "Wheat1", "Right colon"="Goldenrod1","Rectum" = "Sienna1"),
    OS=c("0"="LightYellow", "1" = "Firebrick4"),
    RFS=c("0"="LightYellow", "1" = "Firebrick4"),
    Age = c("white", "firebrick"),
    BMI = c("<25"="LightGoldenrod", "≥25"="Sienna1"),
    Sex = c("Female" = "pink", Male = "LightBlue"),
    DifferentiationGrade = c("Well" = "Wheat1", "Moderate"="Goldenrod1","Poor" = "Sienna1"),
    PathologicalType = c("Mucinous adenocarcinoma"="LightGoldenrod1", "Adenocarcinoma"="Sienna1"),
    HistologicalType = c("Polypoid type" = "LightGoldenrod1", "Ulcerative type" = "Sienna1"),
    TumorSize=c("<5" = "Khaki1", "5-10" = "Goldenrod1",">10" = "Sienna1"),
    TNM.stage=c("Stage I-II" = "Sienna1", "Stage III-IV" = "Firebrick4"),
    Vascular.Thrombosis=c("Negative"= "Azure", "Positive" = "Firebrick4"), 
    PNI=c("Negative"= "Azure", "Positive" = "Firebrick4"), 
    EMVI=c("Negative"= "Azure", "Positive" = "Firebrick4"), 
    Tumor.Budding=c("Negative"= "Azure", "Positive" = "Firebrick4"), 
    BRAF=c("Negative"= "Azure", "Positive" = "Firebrick4"), 
    MSS.MSI=c("dMMR" = "Sienna1", "pMMR" = "Firebrick4")     
)
pheatmap(t(SubtypeInfo),cluster_row=F,clustering_method="ward.D2",show_colnames=F,annotation_col=BeijingSampleInfo, annotation_colors = ann_colors,color = colorRampPalette(c("white","orange","firebrick3"))(50))

pdf("D:/06_CRC/Graph/Part6/BeiJing_clinical.Heatmap.pdf",width=15,height=6)
pheatmap(t(SubtypeInfo),cluster_row=F,clustering_method="ward.D2",show_colnames=F,annotation_col=BeijingSampleInfo, annotation_colors = ann_colors,color = colorRampPalette(c("white","orange","firebrick3"))(50))
dev.off()
pdf("D:/06_CRC/Graph/Part6/BeiJing_clinical.Heatmap_Only.pdf",width=15,height=4)
pheatmap(t(SubtypeInfo),display_numbers =T,number_format = "%.0f",cluster_row=F,clustering_method="ward.D2",show_colnames=F,annotation_colors = ann_colors,color = colorRampPalette(c("white","orange","firebrick3"))(50))
dev.off()

table(BeiJingInfo$GroupByGene)
#NLN_C1 NLN_C2 NLN_C3 NLN_C4 PLN_C1 PLN_C2 
#    50     35     31     32     11     14


t=pheatmap(t(SubtypeInfo),clustering_method="ward.D2",show_colnames=F)
sampleOrder=t$tree_col$label[t$tree_col$order]
SubtypeInfo.df=data.frame(table(BeiJingInfo$MedicalRecordID,BeiJingInfo$CombineGroup))
colnames(SubtypeInfo.df)=c("SampleID","SubType","Number")
SubtypeInfo.df$SampleID=factor(SubtypeInfo.df$SampleID,levels=sampleOrder)
g=ggplot(SubtypeInfo.df,aes(x=SampleID, y=-Number, fill=SubType)) +
      geom_bar(stat = "identity") +
      scale_fill_brewer(palette = "Dark2") +
      #scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_SampleDis_barplot.pdf",width=15,height=2.5)
print(g)
dev.off()


SubTypeList=sort(unique(BeiJingInfo$CombineGroup))
SubTypeNumber=table(BeiJingInfo$CombineGroup)
ShareResult=matrix(data = NA, nrow = length(SubTypeList)*(length(SubTypeList)-1)/2, ncol = 5, byrow = FALSE,dimnames = NULL)
index=1
for(i in 1:(length(SubTypeList)-1)){
      s1=SubTypeList[i]
      s1_Sample=BeiJingInfo[BeiJingInfo$CombineGroup==s1,"MedicalRecordID"]
      for(j in (i+1):length(SubTypeList)){
         s2=SubTypeList[j]
         s2_Sample=BeiJingInfo[BeiJingInfo$CombineGroup==s2,"MedicalRecordID"]
         s1_s2_overlap=intersect(s1_Sample,s2_Sample)
         ShareResult[index,1]=s1
         ShareResult[index,2]=s2
         ShareResult[index,3]=length(s1_s2_overlap)
         ShareResult[index,4]=SubTypeNumber[s1][[1]]
         ShareResult[index,5]=SubTypeNumber[s2][[1]]
         index=index+1
      }
}
ShareResult=data.frame(ShareResult)
colnames(ShareResult)=c("S1","S2","Share","S1N","S2N")
ShareResult$Ratio=as.numeric(ShareResult$Share)/(as.numeric(ShareResult$S1N)+as.numeric(ShareResult$S2N))
ShareResult[order(ShareResult$Ratio,decreasing=T),]
#         S1        S2 Share S1N S2N      Ratio
#1 NLN_C1/C2 NLN_C3/C4    22  85  63 0.14864865
#3 NLN_C3/C4       PLN    11  63  25 0.12500000
#2 NLN_C1/C2       PLN     5  85  25 0.04545455

TargetInfo=c("MedicalRecordID","CombineGroup","LNSize","LNLocation","TumorLocation","Age","Sex","BMI","DifferentiationGrade","PathologicalType","HistologicalType","TumorSize","TNM.stage","MSS.MSI","BRAF","Vascular.Thrombosis","PNI","EMVI","Tumor.Budding","OS","RFS")
BeiJingInfo_SubType=BeiJingInfo[BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2","NLN_C3","NLN_C4"),TargetInfo]
BeiJingInfo_SubType$AgeBin=ifelse(BeiJingInfo_SubType$Age>60,"Old","Young")
BeiJingInfo_SubType$Age=NULL

index=1
p=array()
for(i in 3:length(colnames(BeiJingInfo_SubType))){
   item=colnames(BeiJingInfo_SubType)[i]
   item.mt=as.data.frame.matrix(table(BeiJingInfo_SubType$CombineGroup,BeiJingInfo_SubType[,i]))
   p[index]=fisher.test(item.mt)$p.value
   index=index+1
}
pResult_BeiJing=data.frame(SubType=colnames(BeiJingInfo_SubType)[-c(1:2)],pvalue=p)
pResult_BeiJing=pResult_BeiJing[order(pResult_BeiJing$pvalue),]
g=ggplot(pResult_BeiJing,aes(x=factor(SubType,levels=rev(SubType)), y=-log10(pvalue), fill=-log10(pvalue))) +
      geom_bar(stat = "identity") +
      scale_fill_viridis() +
      theme_bw()+
      coord_flip()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(color="black"),axis.text.y=element_text(color="black"))
pdf("D:/06_CRC/Graph/Part6/BeiJing_CombineGroup_PvalueForAll.barplot.pdf",width=5,height=5)
print(g)
dev.off()
write.table(pResult_BeiJing,file="D:/06_CRC/Graph/Part6/BeiJing_clinical_pvalue.txt",sep="\t",quote=F)



HistologicalType=data.frame(table(BeiJingInfo$CombineGroup,BeiJingInfo$HistologicalType))
colnames(HistologicalType)=c("SubTypes","HistologicalType","Number")
g=ggplot(HistologicalType,aes(x=SubTypes, y=Number, fill=HistologicalType)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("Polypoid type" = "LightGoldenrod1", "Ulcerative type" = "Sienna1")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_CombineGroup_HistologicalType.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
HistologicalType.mt=as.data.frame.matrix(table(BeiJingInfo$CombineGroup,BeiJingInfo$HistologicalType))
fisher.test(HistologicalType.mt[c(1:2),c(1:2)])  #p-value = 0.02004
fisher.test(HistologicalType.mt) #p-value = 0.02784

PNI=data.frame(table(BeiJingInfo$CombineGroup,BeiJingInfo$PNI))
colnames(PNI)=c("SubTypes","PNI","Number")
g=ggplot(PNI,aes(x=SubTypes, y=Number, fill=PNI)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("Azure","Firebrick4")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_CombineGroup_PNI.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
PNI.mt=as.data.frame.matrix(table(BeiJingInfo$CombineGroup,BeiJingInfo$PNI))
fisher.test(PNI.mt[c(1:2),c(1:2)])  #p-value = 0.0387
fisher.test(PNI.mt) #p-value = 0.001209


Vascular.Thrombosis=data.frame(table(BeiJingInfo$CombineGroup,BeiJingInfo$Vascular.Thrombosis))
colnames(Vascular.Thrombosis)=c("SubTypes","Vascular.Thrombosis","Number")
g=ggplot(Vascular.Thrombosis,aes(x=SubTypes, y=Number, fill=Vascular.Thrombosis)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("Azure","Firebrick4")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_CombineGroup_Vascular.Thrombosis.barplot.pdf",width=3.5,height=2.5)
print(g)
dev.off()
Vascular.Thrombosis.mt=as.data.frame.matrix(table(BeiJingInfo$CombineGroup,BeiJingInfo$Vascular.Thrombosis))
fisher.test(Vascular.Thrombosis.mt[c(1:2),c(1:2)])  #p-value = 0.05589
fisher.test(Vascular.Thrombosis.mt) #p-value = 0.009325

TNMStage=data.frame(table(BeiJingInfo$CombineGroup,BeiJingInfo$TNM.stage))
colnames(TNMStage)=c("SubTypes","Stage","Number")
g=ggplot(TNMStage,aes(x=SubTypes, y=Number, fill=Stage)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("Sienna1","Firebrick4")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_CombineGroup_TNMStage.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
TNMStage.mt=as.data.frame.matrix(table(BeiJingInfo$CombineGroup,BeiJingInfo$TNM.stage))
fisher.test(TNMStage.mt[c(1:2),c(1:2)])  #p-value = 0.0406
fisher.test(TNMStage.mt) #p-value = 7.232e-07


OSStage=data.frame(table(BeiJingInfo$CombineGroup,BeiJingInfo$OS))
colnames(OSStage)=c("SubTypes","OSStage","Number")
g=ggplot(OSStage,aes(x=SubTypes, y=Number, fill=OSStage)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("LightYellow","Firebrick4")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_CombineGroup_OSStage.three.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
OSStage.mt=as.data.frame.matrix(table(BeiJingInfo$CombineGroup,BeiJingInfo$OS))
fisher.test(OSStage.mt[c(1:2),c(1:2)])  #p-value = 0.01875
fisher.test(OSStage.mt) #p-value = 0.0001311


RFSStage=data.frame(table(BeiJingInfo$CombineGroup,BeiJingInfo$RFS))
colnames(RFSStage)=c("SubTypes","RFSStage","Number")
g=ggplot(RFSStage,aes(x=SubTypes, y=Number, fill=RFSStage)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("LightYellow","Firebrick4")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_CombineGroup_RFSStage.three.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
RFSStage.mt=as.data.frame.matrix(table(BeiJingInfo$CombineGroup,BeiJingInfo$RFS))
fisher.test(RFSStage.mt[c(1:2),c(1:2)])  #p-value = 0.01875
fisher.test(RFSStage.mt) #p-value = 0.0001311






LNSizeStatus=data.frame(table(BeiJingInfo$CombineGroup,BeiJingInfo$LNSize))
colnames(LNSizeStatus)=c("SubTypes","LNSize","Number")
LNSizeStatus$LNSize=factor(LNSizeStatus$LNSize,levels=c("Small","Large"))
g=ggplot(LNSizeStatus,aes(x=SubTypes, y=Number, fill=LNSize)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("LightBlue","firebrick")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_CombineGroup_LNSize.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
LNSizeStatus.mt=as.data.frame.matrix(table(BeiJingInfo$CombineGroup,BeiJingInfo$LNSize))
fisher.test(LNSizeStatus.mt[c(1:2),c(1:2)])  #p-value = 0.4887
fisher.test(LNSizeStatus.mt) #p-value = 0.06907


LNLocationStatus=data.frame(table(BeiJingInfo$CombineGroup,BeiJingInfo$LNLocation))
colnames(LNLocationStatus)=c("SubTypes","LNLocation","Number")
g=ggplot(LNLocationStatus,aes(x=SubTypes, y=Number, fill=LNLocation)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("LightBlue","orange","Red")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_CombineGroup_LNLocation.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
LNLocationStatus.mt=as.data.frame.matrix(table(BeiJingInfo$CombineGroup,BeiJingInfo$LNLocation))
fisher.test(LNLocationStatus.mt[c(1:3),c(1:3)])  #p-value = 0.1674
fisher.test(LNLocationStatus.mt) #p-value = 0.06907


#The relationship between lymph node size and subtype should be clearly demarcated as its own analysis, rather than being inferred from broader trends.
BeiJingInfo_SubType=BeiJingInfo[BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2","NLN_C3","NLN_C4"),TargetInfo]
GroupByGene_Colors=c("#63B99C","#BB968C","#D986B5","#CAD246","#E2BF8C","#147286")
names(GroupByGene_Colors)=sort(unique(BeiJingInfo$GroupByGene))

LNSizeStatus=data.frame(table(BeiJingInfo$GroupByGene,BeiJingInfo$LNSize))
colnames(LNSizeStatus)=c("SubTypes","LNSize","Number")
LNSizeStatus$LNSize=factor(LNSizeStatus$LNSize,levels=c("Small","Large"))
g=ggplot(LNSizeStatus,aes(x=LNSize, y=Number, fill=SubTypes)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/BeiJing_GroupByGene_LNSize.barplot.pdf",width=3,height=3.5)
print(g)
dev.off()
LNSizeStatus.mt=as.data.frame.matrix(table(BeiJingInfo_SubType$GroupByGene,BeiJingInfo_SubType$LNSize))
fisher.test(LNSizeStatus.mt) #p-value = 0.07643




###########################################################################################
#########################       ShanXi                  ###################################
###########################################################################################

ShanXiInfo=read.table("D:/06_CRC/bulkRNAseq/LNInfo/LymphNodeInfomation-ShanXi_212.txt",header=T,sep="\t")
ShanXiInfo$CombineGroup=ifelse(ShanXiInfo$GroupByGene%in%c("NLN_C1","NLN_C2"),"NLN_C1/C2",ifelse(ShanXiInfo$GroupByGene%in%c("NLN_C3","NLN_C4"),"NLN_C3/C4","PLN"))
ShanXiInfo=ShanXiInfo[ShanXiInfo$GroupByGene%in%c("NLN_C1","NLN_C2","NLN_C3","NLN_C4","PLN_C1"),]
ShanXiSubtypeInfo=as.data.frame.matrix(table(ShanXiInfo$MedicalRecordID,ShanXiInfo$GroupByGene))
length(unique(ShanXiInfo$MedicalRecordID))
pheatmap(t(ShanXiSubtypeInfo),cluster_row=F,clustering_method="ward.D2",color = colorRampPalette(c("white","orange","firebrick3"))(50),show_colnames=F)


TargetInfo=c("MedicalRecordID","Age","Sex","BMI","TumorLocation","PathologicalType","HistologicalType","DifferentiationGrade","TumorSize","TNM.stage","MSS.MSI","BRAF","Vascular.Thrombosis","PNI","EMVI","Tumor.Budding")
ShanXiSampleInfo=unique(ShanXiInfo[,rev(TargetInfo)])
dim(ShanXiSampleInfo)
length(unique(ShanXiSampleInfo$MedicalRecordID))
t=table(ShanXiSampleInfo$MedicalRecordID)
ShanXiSampleInfo[ShanXiSampleInfo$MedicalRecordID%in%names(t[t==2]),]
ShanXiSampleInfo=ShanXiSampleInfo[!is.na(ShanXiSampleInfo$MedicalRecordID),]
rownames(ShanXiSampleInfo)=ShanXiSampleInfo$MedicalRecordID
ShanXiSampleInfo$MedicalRecordID=NULL
# Specify colors
#ShanXiSampleInfo=ShanXiSampleInfo[, apply(ShanXiSampleInfo, 2, function(y) any(!is.na(y)))]
pdf("D:/06_CRC/Graph/Part6/ShanXi_clinical.Heatmap.pdf",width=15,height=6)
pheatmap(t(ShanXiSubtypeInfo),cluster_row=F,clustering_method="ward.D2",show_colnames=F,annotation_col=ShanXiSampleInfo, annotation_colors = ann_colors,color = colorRampPalette(c("white","orange","firebrick3"))(50))
dev.off()

t=pheatmap(t(ShanXiSubtypeInfo),clustering_method="ward.D2",show_colnames=F)
sampleOrder=t$tree_col$label[t$tree_col$order]
ShanXiSubtypeInfo.df=data.frame(table(ShanXiInfo$MedicalRecordID,ShanXiInfo$CombineGroup))
colnames(ShanXiSubtypeInfo.df)=c("SampleID","SubType","Number")
ShanXiSubtypeInfo.df$SampleID=factor(ShanXiSubtypeInfo.df$SampleID,levels=sampleOrder)
g=ggplot(ShanXiSubtypeInfo.df,aes(x=SampleID, y=-Number, fill=SubType)) +
      geom_bar(stat = "identity") +
      scale_fill_brewer(palette = "Dark2") +
      #scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/ShanXi_SampleDis_barplot.pdf",width=15,height=2.5)
print(g)
dev.off()

SubTypeList=sort(unique(ShanXiInfo$CombineGroup))
SubTypeNumber=table(ShanXiInfo$CombineGroup)
ShareResult=matrix(data = NA, nrow = length(SubTypeList)*(length(SubTypeList)-1)/2, ncol = 5, byrow = FALSE,dimnames = NULL)
index=1
for(i in 1:(length(SubTypeList)-1)){
      s1=SubTypeList[i]
      s1_Sample=ShanXiInfo[ShanXiInfo$CombineGroup==s1,"MedicalRecordID"]
      for(j in (i+1):length(SubTypeList)){
         s2=SubTypeList[j]
         s2_Sample=ShanXiInfo[ShanXiInfo$CombineGroup==s2,"MedicalRecordID"]
         s1_s2_overlap=intersect(s1_Sample,s2_Sample)
         ShareResult[index,1]=s1
         ShareResult[index,2]=s2
         ShareResult[index,3]=length(s1_s2_overlap)
         ShareResult[index,4]=SubTypeNumber[s1][[1]]
         ShareResult[index,5]=SubTypeNumber[s2][[1]]
         index=index+1
      }
}
ShareResult=data.frame(ShareResult)
colnames(ShareResult)=c("S1","S2","Share","S1N","S2N")
ShareResult$Ratio=as.numeric(ShareResult$Share)/(as.numeric(ShareResult$S1N)+as.numeric(ShareResult$S2N))
ShareResult[order(ShareResult$Ratio,decreasing=T),]
#         S1        S2 Share S1N S2N      Ratio
#1 NLN_C1/C2 NLN_C3/C4    17 155  40 0.08717949
#3 NLN_C3/C4       PLN     6  40  32 0.08333333
#2 NLN_C1/C2       PLN     9 155  32 0.04812834

HistologicalType=data.frame(table(ShanXiInfo$CombineGroup,ShanXiInfo$HistologicalType))
colnames(HistologicalType)=c("SubTypes","HistologicalType","Number")
HistologicalType=HistologicalType[!HistologicalType$HistologicalType=="",]
g=ggplot(HistologicalType,aes(x=SubTypes, y=Number, fill=HistologicalType)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("Polypoid type" = "LightGoldenrod1", "Ulcerative type" = "Sienna1")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/ShanXi_CombineGroup_HistologicalType.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()

HistologicalType.mt=matrix(HistologicalType[HistologicalType$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4"),"Number"],nrow=2)
fisher.test(HistologicalType.mt)  #p-value = 0.5464
HistologicalType.mt=matrix(HistologicalType[HistologicalType$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4","PLN"),"Number"],nrow=2)
fisher.test(HistologicalType.mt) #p-value =0.2256


PNI=data.frame(table(ShanXiInfo$CombineGroup,ShanXiInfo$PNI))
colnames(PNI)=c("SubTypes","PNI","Number")
PNI=PNI[!PNI$PNI=="",]
g=ggplot(PNI,aes(x=SubTypes, y=Number, fill=PNI)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("Negative"= "Azure", "Positive" = "Firebrick4")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/ShanXi_CombineGroup_PNI.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
PNI.mt=matrix(PNI[PNI$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4"),"Number"],nrow=2)
fisher.test(PNI.mt)  #p-value = 0.001175
PNI.mt=matrix(PNI[PNI$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4","PLN"),"Number"],nrow=3)
fisher.test(PNI.mt) #p-value =0.000999


Vascular.Thrombosis=data.frame(table(ShanXiInfo$CombineGroup,ShanXiInfo$Vascular.Thrombosis ))
colnames(Vascular.Thrombosis)=c("SubTypes","Vascular.Thrombosis","Number")
Vascular.Thrombosis=Vascular.Thrombosis[!Vascular.Thrombosis$Vascular.Thrombosis=="",]
g=ggplot(Vascular.Thrombosis,aes(x=SubTypes, y=Number, fill=Vascular.Thrombosis)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("Negative"= "Azure", "Positive" = "Firebrick4")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/ShanXi_CombineGroup_Vascular.Thrombosis.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
Vascular.Thrombosis.mt=matrix(Vascular.Thrombosis[Vascular.Thrombosis$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4"),"Number"],nrow=2)
fisher.test(Vascular.Thrombosis.mt)  #p-value = 0.03309
Vascular.Thrombosis.mt=matrix(Vascular.Thrombosis[Vascular.Thrombosis$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4","PLN"),"Number"],nrow=3)
fisher.test(Vascular.Thrombosis.mt) #p-value =0.008428

TNMStage=data.frame(table(ShanXiInfo$CombineGroup,ShanXiInfo$TNM.stage))
colnames(TNMStage)=c("SubTypes","Stage","Number")
TNMStage=TNMStage[!TNMStage$Stage=="",]
g=ggplot(TNMStage,aes(x=SubTypes, y=Number, fill=Stage)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("Stage I-II" = "Sienna1", "Stage III-IV" = "Firebrick4")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/ShanXi_CombineGroup_TNMStage.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
TNMStage.mt=matrix(TNMStage[TNMStage$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4"),"Number"],nrow=2)
fisher.test(TNMStage.mt)  #p-value = 0.3155
TNMStage.mt=matrix(TNMStage[TNMStage$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4","PLN"),"Number"],nrow=3)
fisher.test(TNMStage.mt) #4.676e-10


LNSizeStatus=data.frame(table(ShanXiInfo$CombineGroup,ShanXiInfo$LNSize))
colnames(LNSizeStatus)=c("SubTypes","LNSize","Number")
LNSizeStatus=LNSizeStatus[!LNSizeStatus$LNSize=="",]
LNSizeStatus$LNSize=factor(LNSizeStatus$LNSize,levels=c("Small","Large"))
g=ggplot(LNSizeStatus,aes(x=SubTypes, y=Number, fill=LNSize)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("LightBlue","firebrick")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/ShanXi_CombineGroup_LNSizeStatus.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
LNSizeStatus.mt=matrix(LNSizeStatus[LNSizeStatus$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4"),"Number"],nrow=2)
fisher.test(LNSizeStatus.mt)  #p-value = 0.3096
LNSizeStatus.mt=matrix(LNSizeStatus[LNSizeStatus$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4","PLN"),"Number"],nrow=3)
fisher.test(LNSizeStatus.mt) #p-value < 0.1805


LNLocationStatus=data.frame(table(ShanXiInfo$CombineGroup,ShanXiInfo$LNLocation))
colnames(LNLocationStatus)=c("SubTypes","LNLocation","Number")
LNLocationStatus=LNLocationStatus[!LNLocationStatus$LNLocation=="",]
g=ggplot(LNLocationStatus,aes(x=SubTypes, y=Number, fill=LNLocation)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("LightBlue","orange","Red")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/ShanXi_CombineGroup_LNLocation.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
LNLocationStatus.mt=matrix(LNLocationStatus[LNLocationStatus$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4"),"Number"],nrow=2)
fisher.test(LNLocationStatus.mt)  #p-value = 0.3075
LNLocationStatus.mt=matrix(LNLocationStatus[LNLocationStatus$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4","PLN"),"Number"],nrow=3)
fisher.test(LNLocationStatus.mt) #p-value =0.2972

###########################################################################################
#########################       Harbin                  ###################################
###########################################################################################
HarBinInfo=read.table("D:/06_CRC/bulkRNAseq/LNInfo/LymphNodeInfomation-HarBin_162.txt",header=T,sep="\t")
HarBinInfo$CombineGroup=ifelse(HarBinInfo$GroupByGene%in%c("NLN_C1","NLN_C2"),"NLN_C1/C2",ifelse(HarBinInfo$GroupByGene%in%c("NLN_C3","NLN_C4"),"NLN_C3/C4","PLN"))
HarBinInfo=HarBinInfo[HarBinInfo$GroupByGene%in%c("NLN_C1","NLN_C2","NLN_C3","NLN_C4","PLN_C1"),]
HarBinSubtypeInfo=as.data.frame.matrix(table(HarBinInfo$MedicalRecordID,HarBinInfo$GroupByGene))
length(unique(HarBinInfo$MedicalRecordID))
pheatmap(t(HarBinSubtypeInfo),cluster_row=F,clustering_method="ward.D2",color = colorRampPalette(c("white","orange","firebrick3"))(50),show_colnames=F)


ann_colors = list(
    TumorLocation = c("Left colon" = "Wheat1", "Right colon"="Goldenrod1","Rectum" = "Sienna1"),
    Age = c("white", "firebrick"),
    BMI = c("<25"="LightGoldenrod", "≥25"="Sienna1"),
    Sex = c("Female" = "pink", Male = "LightBlue"),
    DifferentiationGrade = c("Well" = "Wheat1", "Moderate"="Goldenrod1","Poor" = "Sienna1"),
    PathologicalType = c("Mucinous adenocarcinoma"="LightGoldenrod1", "Adenocarcinoma"="Sienna1","Squamous cellular carcinoma"="Sienna1"),
    HistologicalType = c("Polypoid type" = "LightGoldenrod1", "Ulcerative type" = "Sienna1"),
    TumorSize=c("<5" = "Khaki1", "5-10" = "Goldenrod1",">10" = "Sienna1"),
    TNM.stage=c("Stage I-II" = "Sienna1", "Stage III-IV" = "Firebrick4"),
    Vascular.Thrombosis=c("Negative"= "Azure", "Positive" = "Firebrick4"), 
    PNI=c("Negative"= "Azure", "Positive" = "Firebrick4"), 
    EMVI=c("Negative"= "Azure", "Positive" = "Firebrick4"), 
    Tumor.Budding=c("Negative"= "Azure", "Positive" = "Firebrick4"), 
    BRAF=c("Negative"= "Azure", "Positive" = "Firebrick4"), 
    MSS.MSI=c("dMMR" = "Sienna1", "pMMR" = "Firebrick4")     
)


TargetInfo=c("MedicalRecordID","Age","Sex","BMI","TumorLocation","PathologicalType","HistologicalType","DifferentiationGrade","TNM.stage","MSS.MSI","BRAF","Vascular.Thrombosis","PNI","EMVI","Tumor.Budding")
HarBinSampleInfo=unique(HarBinInfo[,rev(TargetInfo)])
dim(HarBinSampleInfo)
length(unique(HarBinSampleInfo$MedicalRecordID))
t=table(HarBinSampleInfo$MedicalRecordID)
HarBinSampleInfo[sort(HarBinSampleInfo$MedicalRecordID%in%names(t[t==2])),]

HarBinSampleInfo=HarBinSampleInfo[!is.na(HarBinSampleInfo$MedicalRecordID),]
rownames(HarBinSampleInfo)=HarBinSampleInfo$MedicalRecordID
HarBinSampleInfo$MedicalRecordID=NULL
# Specify colors
HarBinSampleInfo=HarBinSampleInfo[, apply(HarBinSampleInfo, 2, function(y) any(!is.na(y)))]
pdf("D:/06_CRC/Graph/Part6/HarBin_clinical.Heatmap.pdf",width=15,height=5)
pheatmap(t(HarBinSubtypeInfo),cluster_row=F,clustering_method="ward.D2",show_colnames=F,annotation_col=HarBinSampleInfo, annotation_colors = ann_colors,color = colorRampPalette(c("white","orange","firebrick3"))(50))
dev.off()


t=pheatmap(t(HarBinSubtypeInfo),clustering_method="ward.D2",show_colnames=F)
sampleOrder=t$tree_col$label[t$tree_col$order]
HarBinSubtypeInfo.df=data.frame(table(HarBinInfo$MedicalRecordID,HarBinInfo$CombineGroup))
colnames(HarBinSubtypeInfo.df)=c("SampleID","SubType","Number")
HarBinSubtypeInfo.df$SampleID=factor(HarBinSubtypeInfo.df$SampleID,levels=sampleOrder)
g=ggplot(HarBinSubtypeInfo.df,aes(x=SampleID, y=-Number, fill=SubType)) +
      geom_bar(stat = "identity") +
      scale_fill_brewer(palette = "Dark2") +
      #scale_fill_manual(values=GroupByGene_Colors) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/HarBin_SampleDis_barplot.pdf",width=15,height=2.5)
print(g)
dev.off()

SubTypeList=sort(unique(HarBinInfo$CombineGroup))
SubTypeNumber=table(HarBinInfo$CombineGroup)
ShareResult=matrix(data = NA, nrow = length(SubTypeList)*(length(SubTypeList)-1)/2, ncol = 5, byrow = FALSE,dimnames = NULL)
index=1
for(i in 1:(length(SubTypeList)-1)){
      s1=SubTypeList[i]
      s1_Sample=HarBinInfo[HarBinInfo$CombineGroup==s1,"MedicalRecordID"]
      for(j in (i+1):length(SubTypeList)){
         s2=SubTypeList[j]
         s2_Sample=HarBinInfo[HarBinInfo$CombineGroup==s2,"MedicalRecordID"]
         s1_s2_overlap=intersect(s1_Sample,s2_Sample)
         ShareResult[index,1]=s1
         ShareResult[index,2]=s2
         ShareResult[index,3]=length(s1_s2_overlap)
         ShareResult[index,4]=SubTypeNumber[s1][[1]]
         ShareResult[index,5]=SubTypeNumber[s2][[1]]
         index=index+1
      }
}
ShareResult=data.frame(ShareResult)
colnames(ShareResult)=c("S1","S2","Share","S1N","S2N")
ShareResult$Ratio=as.numeric(ShareResult$Share)/(as.numeric(ShareResult$S1N)+as.numeric(ShareResult$S2N))
ShareResult[order(ShareResult$Ratio,decreasing=T),]
#         S1        S2 Share S1N S2N      Ratio
#1 NLN_C1/C2 NLN_C3/C4    19 101  46 0.12925170
#3 NLN_C3/C4       PLN     6  46  15 0.09836066
#2 NLN_C1/C2       PLN     3 101  15 0.02586207

HistologicalType=data.frame(table(HarBinInfo$CombineGroup,HarBinInfo$HistologicalType))
colnames(HistologicalType)=c("SubTypes","HistologicalType","Number")
HistologicalType=HistologicalType[!HistologicalType$HistologicalType=="",]
g=ggplot(HistologicalType,aes(x=SubTypes, y=Number, fill=HistologicalType)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("Polypoid type" = "LightGoldenrod1", "Ulcerative type" = "Sienna1")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/HarBin_CombineGroup_HistologicalType.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()

HistologicalType.mt=matrix(HistologicalType[HistologicalType$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4"),"Number"],nrow=2)
fisher.test(HistologicalType.mt)  #p-value = 0.4682
HistologicalType.mt=matrix(HistologicalType[HistologicalType$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4","PLN"),"Number"],nrow=3)
fisher.test(HistologicalType.mt) #p-value =0.5072


PNI=data.frame(table(HarBinInfo$CombineGroup,HarBinInfo$PNI))
colnames(PNI)=c("SubTypes","PNI","Number")
PNI=PNI[!PNI$PNI=="",]
g=ggplot(PNI,aes(x=SubTypes, y=Number, fill=PNI)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("Negative"= "Azure", "Positive" = "Firebrick4")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/HarBin_CombineGroup_PNI.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
PNI.mt=matrix(PNI[PNI$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4"),"Number"],nrow=2)
fisher.test(PNI.mt)  #p-value = 0.379
PNI.mt=matrix(PNI[PNI$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4","PLN"),"Number"],nrow=3)
fisher.test(PNI.mt) #p-value =0.4777


Vascular.Thrombosis=data.frame(table(HarBinInfo$CombineGroup,HarBinInfo$Vascular.Thrombosis ))
#not include

TNMStage=data.frame(table(HarBinInfo$CombineGroup,HarBinInfo$TNM.stage))
colnames(TNMStage)=c("SubTypes","Stage","Number")
TNMStage=TNMStage[!TNMStage$Stage=="",]
g=ggplot(TNMStage,aes(x=SubTypes, y=Number, fill=Stage)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("Stage I-II" = "Sienna1", "Stage III-IV" = "Firebrick4")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/HarBin_CombineGroup_TNMStage.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
TNMStage.mt=matrix(TNMStage[TNMStage$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4"),"Number"],nrow=2)
fisher.test(TNMStage.mt)  #p-value = 1.236e-05
TNMStage.mt=matrix(TNMStage[TNMStage$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4","PLN"),"Number"],nrow=3)
fisher.test(TNMStage.mt) #4.691e-11

HarBinInfo$LNSizeGroup=ifelse(HarBinInfo$LNSize<1,"Small","Large")
LNSizeStatus=data.frame(table(HarBinInfo$CombineGroup,HarBinInfo$LNSizeGroup))
colnames(LNSizeStatus)=c("SubTypes","LNSize","Number")
LNSizeStatus=LNSizeStatus[!LNSizeStatus$LNSize=="",]
LNSizeStatus$LNSize=factor(LNSizeStatus$LNSize,levels=c("Small","Large"))
g=ggplot(LNSizeStatus,aes(x=SubTypes, y=Number, fill=LNSize)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("LightBlue","firebrick")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/HarBin_CombineGroup_LNSizeStatus.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
LNSizeStatus.mt=matrix(LNSizeStatus[LNSizeStatus$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4"),"Number"],nrow=2)
fisher.test(LNSizeStatus.mt)  #p-value = 0.04409
LNSizeStatus.mt=matrix(LNSizeStatus[LNSizeStatus$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4","PLN"),"Number"],nrow=3)
fisher.test(LNSizeStatus.mt) #p-value=0.02959


LNLocationStatus=data.frame(table(HarBinInfo$CombineGroup,HarBinInfo$LNLocation))
colnames(LNLocationStatus)=c("SubTypes","LNLocation","Number")
LNLocationStatus=LNLocationStatus[!LNLocationStatus$LNLocation=="",]
g=ggplot(LNLocationStatus,aes(x=SubTypes, y=Number, fill=LNLocation)) +
      geom_bar(position = 'fill',stat = "identity") +
      scale_fill_manual(values=c("LightBlue","orange","Red")) +
      theme_bw()+
      theme(panel.spacing = unit(0.05, "cm", data = NULL))+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf("D:/06_CRC/Graph/Part6/HarBin_CombineGroup_LNLocation.barplot.pdf",width=3,height=2.5)
print(g)
dev.off()
LNLocationStatus.mt=matrix(LNLocationStatus[LNLocationStatus$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4"),"Number"],nrow=2)
fisher.test(LNLocationStatus.mt)  #p-value = 0.121
LNLocationStatus.mt=matrix(LNLocationStatus[LNLocationStatus$SubTypes%in%c("NLN_C1/C2","NLN_C3/C4","PLN"),"Number"],nrow=3)
fisher.test(LNLocationStatus.mt) #p-value =0.2647
