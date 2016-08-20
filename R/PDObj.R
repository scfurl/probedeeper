setClass("PDObj", slots=c(
  eset="ExpressionSet",
  ColObj="ColObj",
  LimmaObj="LimmaObj"), package="probedeeper")

# PD<-new("PDObj")
#
# PD@eset<-finaleset.nhp
# PD@ColObj<-ColObj.nhp
# PD@LimmaObj<-LimmaObjCreate(PD@eset, PD@ColObj)
#
# limma.obj<-autoLIMMA(exprs(finaleset.nhp), PD@ColObj@classvec, 1:5, 1:5, pvalue.thresh=0.05, lfc.thresh=1, adjust.method="fdr", method="separate", printdata=FALSE)
# PD@LimmaObj@Contrasts
# MAboxplot("CXCR3", PD)
# debug (MAboxplot)
