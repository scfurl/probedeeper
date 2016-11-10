####PCA.SF####
PD.PCAplot<-function(xdata, ydata, ColObj, title, level=0.95, legendtitle="group"){
#library(ggplot2)
patientcolors<-ColObj@full$line
colors<-ColObj@match$line
classvector<-ColObj@classvec
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
g<-ggplot2::ggplot(data=dfr, ggplot2::aes(x=x, y=y, color=group)) +
  ggplot2::theme_bw()+
  ggplot2::geom_point(size=1.5, alpha=1)+
 # ggplot2::geom_point(data=centroids, ggplot2::aes(x=x, y=y), color=colors, size=3, shape=43) +
  ggplot2::geom_segment(ggplot2::aes(x=x.centroid, y=y.centroid, xend=x, yend=y), color=patientcolors, size=0.2)+
  ggplot2::geom_hline(yintercept = 0)+
  ggplot2::geom_vline(xintercept = 0)+
  ggplot2::theme(axis.ticks.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank())+
  ggplot2::theme(axis.ticks.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank())+
  ggplot2::labs(x=deparse(substitute(xdata)), y=deparse(substitute(ydata)), title=title)+
  ggplot2::theme(legend.key = ggplot2::element_rect(fill = "white"))+
  ggplot2::scale_colour_manual(name=legendtitle, values=colors, guide = ggplot2::guide_legend(override.aes=ggplot2::aes(fill=NA)))
for(i in 1:length(levels(as.factor(dfr_ell$group)))){
  g<-g+ggplot2::geom_path(data=dfr_ell[dfr_ell$group==levels(as.factor(dfr_ell$group))[i],], ggplot2::aes(x=x, y=y), size=0.2, colour=levels(factor(patientcolors, levels = unique(patientcolors)))[i], linetype=1)
  g<-g+ggplot2::geom_polygon(data=dfr_ell[dfr_ell$group==levels(as.factor(dfr_ell$group))[i],], ggplot2::aes(x=x, y=y), fill=levels(factor(patientcolors, levels = unique(patientcolors)))[i], alpha=0.1, linetype=0)}
  g<-g+ ggplot2::geom_point(data=centroids, ggplot2::aes(x=x, y=y), color=colors, size=3, shape=43) +
    ggplot2::theme(legend.position="right")
  #ggplot2::guides(colour=guide_colorbar())
#ggplot2::scale_color_manual("Legend", values=colors)

 return(g)
}
