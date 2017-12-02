
mergeEsets<-function(elist, rma=T){
  if(is.null(names(elist))){names(elist)<-seq(1, length(elist))}
  if(!all(sapply(elist, class) %in% rep("ExpressionSet", length(list)))) stop("Error - not all esets")
  rord<-rownames(exprs(elist[[1]]))
  mats<-lapply(elist, exprs)
  mats<-lapply(mats, "[", rord,)
  #colnames(fData(elist[[2]]))
  #all(as.character(fData(elist[[1]])$ID)==rord)
  fdat<-fData(elist[[1]])
  pdat<-lapply(elist, pData)
  
  colnames_pdat<-lapply(pdat, colnames)
  common<-Reduce(intersect, colnames_pdat)
  pdat_common<-lapply(pdat, "[", , common)
  pdat_final<-do.call("rbind", pdat_common)
  mats_final<-do.call("cbind", mats)
#  mats_final<-rma(mats_final)
  colnames(mats_final)<-make.unique(colnames(mats_final))
  rownames(pdat_final)<-colnames(mats_final)
  
  pdat_final$ID<-unlist(sapply(names(lapply(pdat_common, nrow)), function(x) rep(x, lapply(pdat_common, nrow)[[x]])))
  final_eset<-ExpressionSet(mats_final, phenoData = AnnotatedDataFrame(pdat_final), featureData = AnnotatedDataFrame(fdat))
}