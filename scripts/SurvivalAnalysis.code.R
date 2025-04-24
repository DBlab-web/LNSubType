library(ggplot2)
library(survival)
library(survminer)

BeiJingInfo=read.table("D:/06_CRC/bulkRNAseq/LNInfo/LymphNodeInfomation-Beijing_173.txt",header=T,sep="\t")
BeiJingInfo$CombineGroup=ifelse(BeiJingInfo$GroupByGene%in%c("NLN_C1","NLN_C2"),"NLN_C1/C2",ifelse(BeiJingInfo$GroupByGene%in%c("NLN_C3","NLN_C4"),"NLN_C3/C4","PLN"))
km_trt_fit <- survfit(Surv(Ostime, OS) ~ GroupByGene, data=BeiJingInfo)
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(km_trt_fit,
          pval = TRUE, conf.int = TRUE,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw(), # Change ggplot2 theme
          palette = c("#E7B800", "#2E9FDF","red","blue"))

BeiJingInfoTmp=BeiJingInfo[BeiJingInfo$CombineGroup%in%c("NLN_C1/C2","NLN_C3/C4"),]
km_trt_fit_tmp <- survfit(Surv(Ostime, OS) ~ CombineGroup, data=BeiJingInfoTmp)
# Change color, linetype by strata, risk.table color by strata
g=ggsurvplot(km_trt_fit_tmp,
          pval = TRUE, conf.int = TRUE,
          ylim=c(0.65,1),
          xlim=c(0,40),
          pval.coord = c(0.8, 0.7),
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw(), # Change ggplot2 theme
          palette = c("#249873", "#D05E15"))
pdf("D:/06_CRC/Graph/Part6/BeiJing_CombineGroup_OSStage.survivalplot.pdf",width=6,height=5)
print(g)
dev.off()


km_trt_fit_tmp1 <- survfit(Surv(RFStime, RFS) ~ CombineGroup, data=BeiJingInfoTmp)
# Change color, linetype by strata, risk.table color by strata
g=ggsurvplot(km_trt_fit_tmp1,
          pval = TRUE, conf.int = TRUE,
          risk.table = TRUE, # Add risk table
          ylim=c(0.6,1),
          xlim=c(0,40),
          pval.coord = c(0.8, 0.65),
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw(), # Change ggplot2 theme
          palette = c("#249873", "#D05E15"))
pdf("D:/06_CRC/Graph/Part6/BeiJing_CombineGroup_RFStime.survivalplot.pdf",width=6,height=5)
print(g)
dev.off()