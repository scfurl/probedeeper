###SUBSETTING PD OBJECTS
###THIS FUNCITON SUBSETS THE exprs, ColObj, classvec and reruns the limmaObj providing a complete complement
###of data for a subsetted object

# levels(PD@ColObj@classvec)
# levels<-levels(PD@ColObj@classvec)[c(7,9,12)]
# PDObj<-PD


PD.Subset<-function(PDObj, levels, redoLimma=TRUE, LD=80, element1=NULL,
                    element2=NULL, pvalue.thresh=0.05, lfc.thresh=1, adjust.method="fdr",
                    method="separate", printdata=FALSE){
  if(class(PDObj)!="PDObj"){stop("Input does not appear to be a LimmaObj")}
  if(!all(levels %in% levels(PDObj@ColObj@classvec))){stop("Input incongruent with PDObj Classvec")}
  cv.ind<-which(PDObj@ColObj@classvec %in% levels)
  new.PDObj<-new("PDObj",
                 ColObj=ColObj.Subset(PDObj@ColObj, levels, LD),
                 eset=PDObj@eset[,PDObj@ColObj@classvec %in% levels])
  if(redoLimma==TRUE){
    if(length(levels)<2){stop("Limma Analysis not possible with < 2 contrasts")}
    new.PDObj@LimmaObj<-LimmaObjCreate(new.PDObj@eset, new.PDObj@ColObj, element1, element2, pvalue.thresh, lfc.thresh, adjust.method, method, printdata)
  }
  return(new.PDObj)
}

ColObj.Subset<-function(ColObj, levels, LD=80){
  if(class(ColObj)!="ColObj"){stop("Input does not appear to be a LimmaObj")}
  if(!all(levels %in% levels(PDObj@ColObj@classvec))){stop("Input incongruent with PDObj classvec")}
  cv.ind<-which(ColObj@classvec %in% levels)
  new.cv<-factor(ColObj@classvec[cv.ind])
  new.assign<-ColObj@assign[ColObj@assign$Group %in% levels(new.cv),]
  new.ColObj<-new("ColObj", assign=new.assign, classvec=new.cv)
  return(ColObjInit(new.ColObj, LD))
}

