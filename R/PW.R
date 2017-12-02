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
  ncpus<-as.numeric(detectCores())
  cat(paste0("Running using ", ncpus, " cores"))
  for(i in 1:length(methods)){
    gsaRes[[i]]<-runGSA(geneLevelStats = get(whichstat$stat[which(whichstat$method==methods[i])]),
                        geneSetStat=methods[i],
                        directions=dir,
                        gsc=gsc,
                        nPerm=nPerm,
                        gsSizeLim=gsSizeLim,
                        ncpus=ncpus)
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


extractPianoPW<-function(obj, dir=NULL, rank_threshold = 10){
  if(is.null(dir)){stop("Must pick a direction - distinct_(up or dn), mixed_(up or dn), or non")}
  if(length(dir)>1){
    make_list==TRUE
    if(all(!dir %in% c("distinct_up", "mixed_up", "non", "distinct_dn", "mixed_dn"))){
      stop("Direction entered incorrectly, must use only the following terms:\n
           'distinct_up', 'mixed_up', 'non', 'distinct_dn', 'mixed_dn'")}
    }
  if(!dir %in% c("distinct_up", "mixed_up", "non", "distinct_dn", "mixed_dn")){
    stop("Direction entered incorrectly, must be one or more of the following:\n
         'distinct_up', 'mixed_up', 'non', 'distinct_dn', 'mixed_dn'")}
  listout<-list()
  df<-data.frame()
  #dir<-"distinct_up"
  for(i in 1:length(dir)){
    colI<-which(c("distinct_dn", "mixed_dn", "non", "mixed_dn", "distinct_up") %in% dir[i])
    rowI<-which(ch$rankMat[,colI]<rank_threshold)
    ch$rankMat[rowI,]
    ch$pMat[rowI,]
    ch$nGenesMat[rowI,]
  }
  
  
  ch$rankMat
}