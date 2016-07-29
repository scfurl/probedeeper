####PCA.SF####
PCA.SF<-function(xdata, ydata, classvector, colors, patientcolors,title, level=0.95){
if(length(levels(classvector))!=length(colors)){stop("vector lengths are not matching")}
if(length((classvector))!=length(patientcolors)){stop("vector lengths are not matching")}
#dfr<-as.data.frame(cbind(pca1$li$Axis1, pca1$li$Axis2))
#dfr<-as.data.frame(dfr)
#patientcolors<-patientcolors[order(patientcolors)]
patientcolors<-patientcolors[order(classvector)]
dfr<-as.data.frame(cbind(xdata,ydata))
colnames(dfr)<-c("x", "y")
dfr$group<-as.factor(classvector)
dfr_ell <- data.frame()
centroids <- aggregate(cbind(x,y)~group,data=dfr,mean)
dfr <- merge(dfr,centroids,by="group",suffixes=c("",".centroid"))
for(g in levels(dfr$group)){
  tempdfr<-dfr[dfr$group==g,]
  dfr_ell <- rbind(dfr_ell, cbind(data.frame(with(tempdfr,ellipse::ellipse(as.numeric(cor(x, y)), scale=c(sd(x),sd(y)),
                                                               centre=c(mean(x),mean(y)), level=level))),group=g))
}
g<-ggplot2::ggplot(data=dfr, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(colour=patientcolors, size=1.5, alpha=1)+
  ggplot2::geom_point(data=centroids, ggplot2::aes(x=x, y=y), colour=colors, size=3, shape=43) +
  ggplot2::geom_segment(ggplot2::aes(x=x.centroid, y=y.centroid, xend=x, yend=y), colour=patientcolors, size=0.2)+
  ggplot2::geom_hline(yintercept = 0)+
  ggplot2::geom_vline(xintercept = 0)+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  ggplot2::theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  ggplot2::labs(x=deparse(substitue(xdata)), y=deparse(substitue(ydata)), title=title)
for(i in 1:length(levels(as.factor(dfr_ell$group)))){
  g<-g+ggplot2::geom_path(data=dfr_ell[dfr_ell$group==levels(as.factor(dfr_ell$group))[i],], ggplot2::aes(x=x, y=y), size=0.2, colour=levels(factor(patientcolors, levels = unique(patientcolors)))[i], linetype=1)
  g<-g+ggplot2::geom_polygon(data=dfr_ell[dfr_ell$group==levels(as.factor(dfr_ell$group))[i],], ggplot2::aes(x=x, y=y), fill=levels(factor(patientcolors, levels = unique(patientcolors)))[i], alpha=0.1, linetype=0) }

print(g)
return(g)
}
