###Integrating piano PW analysis

#PDObj<-PD
#PD@LimmaObj@Contrasts$meaning
#comparison="KY1005.Siro.FvHealthy.Control"
PreConsensusHM<-function(PDObj, gsc, comparison,
                         methods=c("mean", "median", "sum", "maxmean","fisher", "stouffer","tailStrength"),
                         gsSizeLim=c(10,800), nPerm=1000){
  index<-which(PDObj@LimmaObj@Contrasts$meaning %in% comparison)
  inputP<-PDObj@LimmaObj@DEGenes[[index]][,][,"adj.P.Val", drop=FALSE]
  inputT<-PDObj@LimmaObj@DEGenes[[index]][,][,"t", drop=FALSE]
  dir<-PDObj@LimmaObj@DEGenes[[index]][,][,"logFC", drop=FALSE]

  whichstat<-data.frame(method=c("mean", "median", "sum", "maxmean","fisher", "stouffer","tailStrength"),
                        stat=c("inputT", "inputT", "inputT", "inputT", "inputP", "inputP", "inputP"), stringsAsFactors=F)
  gsaRes<-list()
  ncpus<-detectCores()
  for(i in 1:length(methods)){
    gsaRes[[i]]<-runGSA(geneLevelStats = get(whichstat$stat[which(whichstat$method==methods[i])]),
                        geneSetStat=methods[i],
                        directions=dir,
                        gsc=gsc,
                        nPerm=nPerm,
                        gsSizeLim=gsSizeLim,
                        ncpus=4)
  }
  names(gsaRes)<-methods
  return(gsaRes)
}


PDrunGSA<-function(PDObj, gsc, comparison, geneSetStat="fisher", gsSizeLim=c(5,800), nPerm=1000){
  index<-which(PDObj@LimmaObj@Contrasts$meaning %in% comparison)
  inputP<-PDObj@LimmaObj@DEGenes[[index]][,][,"adj.P.Val", drop=FALSE]
  inputT<-PDObj@LimmaObj@DEGenes[[index]][,][,"t", drop=FALSE]
  dir<-PDObj@LimmaObj@DEGenes[[index]][,][,"logFC", drop=FALSE]
  whichstat<-data.frame(method=c("mean", "median", "sum", "maxmean","fisher", "stouffer","tailStrength"),
                        stat=c("inputT", "inputT", "inputT", "inputT", "inputP", "inputP", "inputP"), stringsAsFactors=F)
  Res<-runGSA(geneLevelStats = get(whichstat$stat[which(whichstat$method==geneSetStat)]),
                        geneSetStat=geneSetStat,
                        directions=dir,
                        gsc=gsc,
                        nPerm=nPerm,
                        gsSizeLim=gsSizeLim,
                        ncpus=4)
  return(Res)
}
