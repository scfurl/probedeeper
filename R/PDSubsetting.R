###SUBSETTING PD OBJECTS
###THIS FUNCITON SUBSETS THE exprs, ColObj, classvec and reruns the limmaObj providing a complete complement
###of data for a subsetted object

# levels(PD@ColObj@classvec)
# levels<-levels(PD@ColObj@classvec)[c(7,9,12)]
# PDObj<-PD


PD.Subset<-function(PDObj, bylevels=NULL, bycolumns=NULL, redoLimma=TRUE, LD=80, element1=NULL,
                    element2=NULL, pvalue.thresh=0.05, lfc.thresh=1, adjust.method="fdr",
                    method="separate", printdata=FALSE, pie=F){
  if(class(PDObj)!="PDObj"){stop("Input does not appear to be a LimmaObj")}
  if(!all(bylevels %in% levels(PDObj@ColObj@classvec))){stop("Input incongruent with PDObj Classvec")}
  if(is.null(bylevels) && is.null(bycolumns)){stop("No subsetting paramater set, use bylevels or bycolumns")}
  if(!is.null(bylevels) && !is.null(bycolumns)){stop("Cannot use both bylevels and bycolumns, choose one")}
  if(!is.null(bylevels)){
    cv.ind<-which(PDObj@ColObj@classvec %in% bylevels)
    new.PDObj<-new("PDObj",
                 ColObj=ColObj.Subset(PDObj@ColObj, bylevels, LD, pie=pie),
                 eset=PDObj@eset[,PDObj@ColObj@classvec %in% bylevels])
    if(redoLimma==TRUE){
      if(length(bylevels)<2){stop("Limma Analysis not possible with < 2 contrasts")}
      new.PDObj@LimmaObj<-LimmaObjCreate(new.PDObj@eset, new.PDObj@ColObj, element1, element2, pvalue.thresh, lfc.thresh, adjust.method, method, printdata)
    }
  }
  if(!is.null(bycolumns)){

    new.PDObj<-new("PDObj",
                   ColObj=ColObj.Subset(PDObj@ColObj, levels, LD, pie=pie),
                   eset=PDObj@eset[,PDObj@ColObj@classvec %in% bylevels])
    if(redoLimma==TRUE){
      if(length(levels)<2){stop("Limma Analysis not possible with < 2 contrasts")}
      new.PDObj@LimmaObj<-LimmaObjCreate(new.PDObj@eset, new.PDObj@ColObj, element1, element2, pvalue.thresh, lfc.thresh, adjust.method, method, printdata)
    }
  }
  return(new.PDObj)
}

debug(ColObj.Subset)


ColObj.Subset<-function(ColObj, levels, LD=80, pie=F){
  if(class(ColObj)!="ColObj"){stop("Input does not appear to be a LimmaObj")}
  if(!all(levels %in% levels(ColObj@classvec))){stop("Input incongruent with PDObj classvec")}
  cv.ind<-which(ColObj@classvec %in% levels)
  new.cv<-factor(ColObj@classvec[cv.ind])
  new.assign<-ColObj@assign[ColObj@assign$Group %in% levels(new.cv),]
  new.ColObj<-new("ColObj", assign=new.assign, classvec=new.cv)
  if(pie==T){
    return(ColObjInit(new.ColObj, LD))
  }
  if(pie==F){
    return(ColObjInit(new.ColObj, LD, pie=F))
  }
}

ColObj.Subsetbycols<-function(ColObj, indices, LD=80, pie=F){
  if(class(ColObj)!="ColObj"){stop("Input does not appear to be a LimmaObj")}
  cv.ind<-indices
  new.cv<-factor(as.character(ColObj@classvec[cv.ind]))
  new.assign<-ColObj@assign[ColObj@assign$Group %in% levels(new.cv),]
  new.ColObj<-new("ColObj", assign=new.assign, classvec=new.cv)
  if(pie==T){
    return(ColObjInit(new.ColObj, LD))
  }
  if(pie==F){
    return(ColObjInit(new.ColObj, LD, pie=F))
  }
}



WriteGPFiles<-function(PD, directory, prefix){
  ##This function writes GCT and CLS files for a PD object
  alldata<-exprs(PD@eset)
  alldata.gct<-alldata[,order(colnames(alldata))]
  gsea.write.gct(alldata.gct, file.path(directory, paste(prefix, ".gct")))
  classvec.cls<-PD@ColObj@classvec[order(colnames(alldata))]
  gsea.write.cls(classvec.cls, file.path(directory, paste(prefix, ".cls")))
}

