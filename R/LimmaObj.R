
setClass("LimmaObj", slots=c(
  Contrasts="data.frame",
  DEGenes="list",
  AllGenes="list",
  Inputs="list"), package="probedeeper")


LimmaObjCreate<-function(eset, ColObj, element1=NULL, element2=NULL, pvalue.thresh=0.05, lfc.thresh=1, adjust.method="fdr", method="separate", printdata=FALSE){
  ####Auto-LIMMA####
  #inputs data.sel, classvec.sel, element1,(2), pvalue.thresh, lfc.thresh
  classvec.sel<-ColObj@classvec
  if(is.null(element1) && is.null(element2)){
    element1<-1:length(levels(classvec.sel))
    element2<-element1}
  data.sel<-exprs(eset)
  limmaDE.output<-list()
  mydesign <- model.matrix(~0 + classvec.sel)
  LETTERS2600 <- c(sapply(LETTERS, function(x) paste0(x, 1:100)))
  colnames(mydesign) <- LETTERS2600[seq(from=1, to=length(levels(classvec.sel)))]
  names(colnames(mydesign))<-levels(classvec.sel)
  fit <- limma::lmFit(data.sel, mydesign)
  mydesign<<-mydesign
  contrast.string<-contrast.eq(colnames(mydesign), element1=element1, element2=element2)
  contrast.matrix<-do.call(makeContrasts.SF, as.list(contrast.string))
  space1<-unlist(lapply(gregexpr(" ", contrast.string), "[[", 1))
  space2<-unlist(lapply(gregexpr(" ", contrast.string), "[[", 2))
  space3<-unlist(lapply(gregexpr(" ", contrast.string), "[[", 3))
  space4<-unlist(lapply(gregexpr(" ", contrast.string), "[[", 4))
  colnames(contrast.matrix)<-substring(contrast.string,1,space1-1)
  genofits <- limma::contrasts.fit(fit, contrast.matrix)
  geno_ebFit <- limma::eBayes(genofits)
  results<-limma::decideTests(geno_ebFit, method)
  contrast.list<-as.list(substring(contrast.string,1,space1-1))
  #contrast.list<-as.list(substring(contrast.string,1,2))
  names(contrast.list)<-as.list(substring(contrast.string,1,space1-1))
  contrast.mat<-as.data.frame(mapply(append, contrast.list, substring(contrast.string,space2+1, space3-1)), stringsAsFactors=FALSE)
  contrast.mat<-rbind(contrast.mat, substring(contrast.string,space4+1))
  colnames(contrast.mat)<-paste("Contrast", as.character(1:ncol(contrast.mat)), sep="")
  tmp<-paste(names(colnames(mydesign))[match(as.character(contrast.mat[2,]), colnames(mydesign))],
             names(colnames(mydesign))[match(contrast.mat[3,], colnames(mydesign))], sep="v")
  contrast.mat<-rbind(contrast.mat, tmp)
  # contrast.list<-list(coef=contrast.mat[1,], meaning=contrast.mat[4,])
  contrast.df<-data.frame(label=colnames(contrast.mat), coef=as.character(contrast.mat[1,]), meaning=as.character(contrast.mat[4,]))
  #limmaDE.output<-contrast.list
  LimmaObj<-new('LimmaObj')
  LimmaObj@Contrasts<-contrast.df
  tmp.list<-list()
  tmp.list2<-list()
  for(i in 1:length(contrast.list[[1]])){
    tmp.list2[[i]]<-limma::topTable(geno_ebFit,coef=as.character(limmaDE.output[[1]][i]),n=60000,adjust.method=adjust.method, sort.by="logFC", p.value=pvalue.thresh,lfc=lfc.thresh)
    names(tmp.list2)[i]<-paste(as.character(limmaDE.output[[1]][i]), "-", as.character(limmaDE.output[[2]][i]), sep="")
    #if(os=="Windows"){rownames(tmp.list2[[i]])<-tmp.list2[[i]]$ID}
    tmp.list[[i]]<-limma::topTable(geno_ebFit,coef=as.character(limmaDE.output[[1]][i]),n=60000,adjust.method=adjust.method, sort.by="logFC")
    names(tmp.list)[i]<-paste(as.character(limmaDE.output[[1]][i]), "-", as.character(limmaDE.output[[2]][i]), sep="")
    #if(os=="Windows"){rownames(tmp.list[[i]])<-tmp.list[[i]]$ID}
  }
  LimmaObj@DEGenes<-tmp.list2
  LimmaObj@AllGenes<-tmp.list
#   names(limmaDE.output)[4]<-"AllLimmaDEData"
#   names(limmaDE.output)[3]<-"SigLimmaDEData"
#   comp1<-names(colnames(mydesign))[match(contrast.mat[2,], colnames(mydesign))]
#   comp2<-names(colnames(mydesign))[match(contrast.mat[3,], colnames(mydesign))]
#   limmaDE.output[5]<-list(Comparisons.df=data.frame(Comp1=comp1, Comp2=comp2))
  LimmaObj@Inputs<-list(p.adjust=adjust.method, p.thresh=pvalue.thresh, lfc.thresh=lfc.thresh)
#   if(printdata==TRUE){
#     print(names(limmaDE.output[[3]]))
#   }
  return(LimmaObj)
}
