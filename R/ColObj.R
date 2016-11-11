##ColObj##
#This set of functions allows for easy color reference for microarray analysis.
#Essentially, the ColObj has three slots 1) an assignment slot 2) a slot for full color
#reference (i.e for every sample) and 3) a match slot for an ordered list of colors by matching classvec

ColObjInit<-function(ColObj, LD=80, pie-T){
  if(class(ColObj)!="ColObj"){stop("Input is not a ColObj class")}
  selected<-levels(ColObj@classvec)
  ColObj@assign<-ColObj@assign[match(selected,ColObj@assign$Group),]
  colmatch<-vector()
  if(substring(ColObj@assign$Color[1], 1,1)!="#") {
    colmatch<-paste("#", ColObj@assign$Color[match(selected,ColObj@assign$Group)], sep="")
  }
  else
  {
    colmatch<-ColObj@assign$Color[match(selected,ColObj@assign$Group)]
  }
  names(colmatch)<-ColObj@assign$Group[match(selected,ColObj@assign$Group)]
  ColObj@full$line<-colmatch[match(ColObj@classvec, names(colmatch))]
  ColObj@full$fill<-sapply(ColObj@full$line, LightenDarkenColor, LD)
  ColObj@match$line<-colmatch
  ColObj@match$fill<-sapply(colmatch, LightenDarkenColor, LD)
  if(pie=T){
  par(mfrow=c(1,2))
  pie(rep(1,length(ColObj@match$line)), col=ColObj@match$line, labels=names(ColObj@match$line), main="Line Colors")
  pie(rep(1,length(ColObj@match$fill)), col=ColObj@match$fill, labels=names(ColObj@match$fill), main="Fill Colors")
  par(mfrow=c(1,1))
  }
  return(ColObj)
}



ColObjCreate<-function(classvec, LD=80, type='rainbow', pie=TRUE){
  if(type=='rainbow'){
    if(is.null(classvec)==FALSE & class(classvec)=="factor"){
#     file<-"~/Dropbox (Kean Lab)/AWS/Scott/Rproj/Kymab/Classvec"
#     classvec<-readRDS(file)
    selected<-levels(classvec)
    colmatch<-substring(rainbow(length(selected)), 1, 7)
    names(colmatch)<-selected
    color.ws<-data.frame(Group=selected, ColorLabel=make.unique(rep("Color", length(selected))), Color=colmatch)
    ColObj<-new('ColObj', assign=color.ws, classvec=classvec)
    ColObj@full$line<-colmatch[match(ColObj@classvec, names(colmatch))]
    ColObj@full$fill<-sapply(ColObj@full$line, LightenDarkenColor, LD)
    ColObj@match$line<-colmatch
    ColObj@match$fill<-sapply(colmatch, LightenDarkenColor, LD)
    if(pie==TRUE){
    par(mfrow=c(1,2))
    pie(rep(1,length(ColObj@match$line)), col=ColObj@match$line, labels=names(ColObj@match$line), main="Line Colors")
    pie(rep(1,length(ColObj@match$fill)), col=ColObj@match$fill, labels=names(ColObj@match$fill), main="Fill Colors")
    par(mfrow=c(1,1))
    }
    return(ColObj)
    }
    else{stop("Check Input - not correct - either not a factor or null")}
  }
}

LightenDarkenColor<-function(col, amt) {
  if (substring(col, 1, 1)=="#") {
    col = substring(col, 2)
  }
  num = as.hexmode(col)
  r = bitwShiftR(num, 16) + amt
  if (r > 255) {r = 255}
  if  (r < 0) {r = 0}
  b = bitwAnd(bitwShiftR(num, 8), 0x00FF) + amt
  if (b > 255) {b = 255}
  if  (b < 0) {b = 0}
  g = bitwAnd(num, 0x0000FF) + amt
  if (g > 255) {g = 255}
  if (g < 0) {g = 0}
  inter<-paste("000000", as.hexmode(bitwOr(g , bitwOr(bitwShiftL(b, 8), bitwShiftL(r, 16)))), sep="")
  ret<-substr(inter, nchar(inter)-5, nchar(inter))
  return(toupper(paste("#", ret, sep="")))
}
