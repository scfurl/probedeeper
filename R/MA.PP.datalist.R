###analyseMAdatalist###



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


