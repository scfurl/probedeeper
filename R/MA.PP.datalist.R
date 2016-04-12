###analyseMAdatalist###
setClass("MA.PP.datalist", representation(
  Suffix="character",
  Phenodata="data.frame",
  Data="list",
  Other="list"))


analyseMAdatalist<-function(data.list){
  pdf(file=spaste(format(Sys.time(), "%y%m%d-"), data.list@Suffix, "-output.pdf"), height=8.5, width=11)
  for(i in 1:length(data.list@Data)){
    par(mfrow = c(2,2), cex=0.5, xpd=TRUE)
    data.tmp<-data.list@Data[[i]]
    gcolors<-data.list@Other$gcolors
    title<-names(data.list@Data)[i]
    SFHist(data.tmp, gcolors, paste("Histogram of ", title, "\n", data.list@Suffix, sep=""))
    boxplot(data.tmp, col=gcolors, main=paste("Histogram of ", title, "\n", data.list@Suffix, sep=""))
    par(cex=0.75)
    colnames(data.tmp)<-make.unique(colnames(data.tmp))
    pca1 <- dudi.pca(t(data.tmp), scann = FALSE, nf=10)
    s.class(dfxy = pca1$li, fac = data.list@Other$classvec, col=data.list@Other$gcolor, xax = 1, yax = 2,
            sub=paste("Groups-", data.list@Suffix, names(data.list@Data)[i], sep="-"))
    s.class(dfxy = pca1$li, fac = as.factor(data.list@Phenodata$Batch), col=c("darkgreen", "darkorange"), xax = 1, yax = 2, 
            sub=paste("Batches-", data.list@Suffix, names(data.list@Data)[i], sep="-"))
  }
  dev.off()
}

SFHist<-function(matrix, colors=rainbow(ncol(matrix)), title){
  density.data<-apply(matrix, 2, density)
  yrange<-max(unlist(lapply(density.data, function(x) max(x$y))))
  ncols=ncol(matrix)
  plot(density(matrix[,1]), col=colors[1], main=title, ylim=c(0,yrange))
  for(i in 2:ncols){
    lines(density(matrix[,i]), col=colors[i])
  }
}

getColors<-function(classvec, colFile, pie=FALSE){
  if(file.exists(colFile)==FALSE){stop("file does not exist")}
  pal<-read.table(colFile, header=TRUE, stringsAsFactors=FALSE)
  cols<-list(cols=c(rgb(pal$R/255,pal$G/255,pal$B/255))[!is.na(pal$Group)], name=pal$Name[!is.na(pal$Group)], group=pal$Group[!is.na(pal$Group)])
  selected<-levels(as.factor(classvec))
  colmatch<-cols$cols[match(selected,cols$group)]
  names(colmatch)<-cols$group[match(selected,cols$group)]
  #pie(rep(1,length(colmatch)), col=colmatch, labels=names(colmatch))
  if(pie==TRUE){
    pie(rep(1,length(gcolor)), col=gcolor, labels=names(gcolor))
  }
  gcolor<-colmatch[order(names(colmatch))]
  return(gcolor)
}
