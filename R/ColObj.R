##ColObj##
#This set of functions allows for easy color reference for microarray analysis.
#Essentially, the ColObj has three slots 1) an assignment slot 2) a slot for full color
#reference (i.e for every sample) and 3) a match slot for an ordered list of colors by matching classvec

setClass("ColObj", slots=c(
  assign="data.frame",
  classvec="factor",
  full="list",
  match="list"), package="probedeeper")


ColObjInit<-function(ColObj, LD=80){
  if(class(ColObj)!="ColObj"){stop("Input is not a ColObj class")}
  selected<-levels(ColObj@classvec)
  colmatch<-paste("#", ColObj@assign$Color[match(selected,ColObj@assign$Group)], sep="")
  names(colmatch)<-ColObj@assign$Group[match(selected,ColObj@assign$Group)]
  ColObj@full$line<-colmatch[match(ColObj@classvec, names(colmatch))]
  ColObj@full$fill<-sapply(ColObj@full$line, LightenDarkenColor, LD)
  ColObj@match$line<-colmatch
  ColObj@match$fill<-sapply(colmatch, LightenDarkenColor, LD)
  par(mfrow=c(1,2))
  pie(rep(1,length(ColObj@match$line)), col=ColObj@match$line, labels=names(ColObj@match$line), main="Line Colors")
  pie(rep(1,length(ColObj@match$fill)), col=ColObj@match$fill, labels=names(ColObj@match$fill), main="Fill Colors")
  return(ColObj)
}

ColObjCreate<-function(classvec, LD=80, type='rainbow'){
  if(type=='rainbow'){
#     file<-"~/Dropbox (Kean Lab)/AWS/Scott/Rproj/Kymab/Classvec"
#     classvec<-readRDS(file)
    selected<-levels(classvec)
    colmatch<-rainbow(length(selected))
    names(colmatch)<-selected
    color.ws<-data.frame(Group=selected, ColorLabel=make.unique(rep("Color", length(selected))), Color=    sapply(colmatch, substring, 2,9))
    ColObj<-new('ColObj', assign=color.ws, classvec=classvec)
    ColObj@full$line<-colmatch[match(ColObj@classvec, names(colmatch))]
    ColObj@full$fill<-sapply(ColObj@full$line, LightenDarkenColor, LD)
    ColObj@match$line<-colmatch
    ColObj@match$fill<-sapply(colmatch, LightenDarkenColor, LD)
    par(mfrow=c(1,2))
    pie(rep(1,length(ColObj@match$line)), col=ColObj@match$line, labels=names(ColObj@match$line), main="Line Colors")
    pie(rep(1,length(ColObj@match$fill)), col=ColObj@match$fill, labels=names(ColObj@match$fill), main="Fill Colors")
  }


}

