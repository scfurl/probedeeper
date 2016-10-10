MAboxplot<-function(gene=NULL,
                    PDObj,
                    alpha=0.8,
                    dot.size=6,
                    box.size=1,
                    box.width=1,
                    reorder=NULL,
                    stat.test="limma",
                    annotate=TRUE,
                    p.value=0.05,
                    sampleNames=NULL,
                    error.bars=TRUE){
  if(class(PDObj)!="PDObj"){stop("PDObj does not appear correct")}
  if(is.null(gene)==TRUE){stop("No gene input")}
  ColObj<-PDObj@ColObj
  classvec<-as.factor(ColObj@classvec)
  line.cols<-ColObj@match$line
  dot.fill.cols<-ColObj@match$fill
  box.fill.cols<-ColObj@match$fill
  if(is.null(reorder)==FALSE){classvec<-factor(classvec,levels(classvec)[reorder])
  line.cols<-line.cols[reorder]
  dot.fill.cols<-dot.fill.cols[reorder]
  box.fill.cols<-box.fill.cols[reorder]
  }
  if(is.null(ColObj)==TRUE){classvec<-factor(classvec,levels(classvec)[reorder])
  line.cols<-line.cols[reorder]
  dot.fill.cols<-dot.fill.cols[reorder]
  box.fill.cols<-box.fill.cols[reorder]
  }
  eset<-PDObj@eset
  obj<-Biobase::exprs(eset)[gene,]
  df<-data.frame(gp=classvec, y=obj)
  maxh <- max(df$y)
  minh<-min(df$y)
  spread<-(maxh-minh)/14
  maxh<-maxh+spread
  df$maxh<-maxh
  df$spread<-spread
  ds <- plyr::ddply(df, plyr::as.quoted("gp"), plyr::summarise, mean = mean(y), sd = sd(y))
  if(stat.test=="pairwiset"){
    stats<-pairwise.t.test(df$y, classvec, p.adjust="fdr")
  }
  if(stat.test=="limma"){

    stat.df<-data.frame(Comp1=as.factor(PDObj@LimmaObj@Contrasts$comp1), Comp2=as.factor(PDObj@LimmaObj@Contrasts$comp2), p.value=ExtractLimmaObj(gene, PDObj@LimmaObj)$adj.P.Val, stringsAsFactors=FALSE)
    stat.df<-stats.table(stat.df)
    if(length(levels(classvec))==2){
      if(is.null(reorder)==FALSE){
        stat.df<-t(data.frame(stat.df[reorder,]))
        colnames(stat.df)<-levels(classvec)[reorder]
      }
      else{
        stat.df<-t(data.frame(stat.df))
      }
      rownames(stat.df)<-levels(classvec)[2]
    }
    else{
      if(is.null(reorder)==FALSE){stat.df<-stat.df[,reorder]
      }
      rownames(stat.df)<-levels(classvec)[order(levels(classvec))][2:length(levels(classvec))]
    }

    #stat.df[is.na(stat.df)]<-1
    stats<-list(p.value=stat.df)
    array.ind<-as.data.frame(which(stats$p.value < p.value, arr.ind=T))
    if(length(levels(classvec))==2){
      if(is.null(reorder)==TRUE){
        array.ind$row<-match(rownames(stat.df)[array.ind$row], levels(classvec))
      }
    }
    else{array.ind$row<-match(rownames(stat.df)[array.ind$row], levels(classvec))}
    rownames(array.ind)<-NULL
    colnames(array.ind)<-c("start","end")
    array.ind$y<-seq(from=maxh+0.5*(spread), by=spread/2, length.out=nrow(array.ind))
  }
  if(stat.test=="pairwiset"){
    array.ind<-as.data.frame(which(stats$p.value < p.value, arr.ind=T))
    rownames(array.ind)<-NULL
    array.ind$row<-array.ind$row+1
    colnames(array.ind)<-c("start","end")
    array.ind$y<-seq(from=maxh+0.5*(spread), by=spread/2, length.out=nrow(array.ind))
  }
  #add sampleNames to df
  if(length(sampleNames)!=0){
    df$sampleNames<-sampleNames
    g<-ggplot2::ggplot(df, ggplot2::aes(x = gp, y = y)) +
      ggplot2::geom_boxplot(size=box.size, alpha=0.6, fill=box.fill.cols, colour=line.cols, width=box.width, outlier.size=NULL, width=box.width) +
      ggplot2::geom_point(size=dot.size, shape=21, colour=line.cols[classvec], fill=dot.fill.cols[classvec], alpha=alpha, position = ggplot2::position_jitter(width = .1)) +
      ggplot2::geom_text(data=df, ggplot2::aes(x = gp, y = y, label=sampleNames), size = 3, hjust=-1) +
      ggplot2::labs(list(x = NULL, y = "Log2 Transformed Data", title=gene)) +
      #ylab(expression(paste("Log", [2], " Transformed Data", sep="")))+
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text=ggplot2::element_text(size=16),
                     axis.title.x=ggplot2::element_text(size=20, vjust=1),
                     axis.text.x = ggplot2::element_text(angle = 45, hjust=1),
                     axis.title.y=ggplot2::element_text(size=16, vjust=0.5), plot.title = ggplot2::element_text(vjust = 0, size=20),
                     axis.line = ggplot2::element_line(colour = "black"),
                     #text=ggplot2::element_text(family="Myriad Pro"),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
    if(nrow(array.ind)==0){
      print(g)
    }
    else
    {if(error.bars==TRUE){
      for(i in 1:nrow(array.ind)){
        g<-g+ggplot2::geom_segment(ggplot2::aes_string(x = array.ind$start[i], y = array.ind$y[i], xend = array.ind$end[i], yend=array.ind$y[i]), lwd=0.5,arrow = arrow(angle=90, ends="both", length = grid::unit(0.1, "cm")))
      }
    }
    }
    print(g)
    footie1<-ifelse(stat.test=="limma", paste("Mod. Bayesian T statistic corrected using ", PD@LimmaObj@Inputs$p.adjust, sep=""), "Pairwise T Test, FDR-corrected")
    Footnote.txt<-paste("Horizontal bars indicate p <0.05 using ", footie1, sep="")
    makeFootnote(Footnote.txt,  color = "black")
  }
  else
  {
    g<-ggplot2::ggplot(df, ggplot2::aes(x = gp, y = y)) +
      ggplot2::geom_boxplot(size=box.size, alpha=0.6, fill=box.fill.cols, width=box.width, colour=line.cols, outlier.size=NULL) +
      ggplot2::geom_point(size=dot.size, shape=21, colour=line.cols[classvec], fill=dot.fill.cols[classvec], alpha=alpha, position = ggplot2::position_jitter(width = .1)) +
      ggplot2::labs(list(x = NULL, y = "Log2 Transformed Data", title=gene)) +
      #ylab(expression(paste("Log", [2], " Transformed Data", sep="")))+
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text=ggplot2::element_text(size=16),
                     axis.title.x=ggplot2::element_text(size=20, vjust=1),
                     axis.text.x = ggplot2::element_text(angle = 45, hjust=1),
                     axis.title.y=ggplot2::element_text(size=16, vjust=0.5), plot.title = ggplot2::element_text(vjust = 0, size=20),
                     axis.line = ggplot2::element_line(colour = "black"),
                     #text=ggplot2::element_text(family="Myriad Pro"),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
    if(nrow(array.ind)==0){
    }
    else
    {
      for(i in 1:nrow(array.ind)){
        g<-g+ggplot2::geom_segment(ggplot2::aes_string(x = array.ind$start[i], y = array.ind$y[i], xend = array.ind$end[i], yend=array.ind$y[i]), lwd=0.5,arrow = ggplot2::arrow(angle=90, ends="both", length = grid::unit(0.1, "cm")))
      }
    }
    print(g)
    footie1<-ifelse(stat.test=="limma", paste("Mod. Bayesian T statistic corrected using ", PDObj@LimmaObj@Inputs$p.adjust, sep=""), "Pairwise T Test, FDR-corrected")
    Footnote.txt<-paste("Horizontal bars indicate p <0.05 using ", footie1, sep="")
    makeFootnote(Footnote.txt,  color = "black")
  }
}

# write a simple function to add footnote
makeFootnote <- function(footnoteText =
                           format(Sys.time(), "%d %b %Y"),
                         size = .7, color = grey(.5))
{
  grid::pushViewport(grid::viewport())
  grid::grid.text(label = footnoteText ,
                  x = grid::unit(1,"npc") - grid::unit(2, "mm"),
                  y = grid::unit(2, "mm"),
                  just = c("right", "bottom"),
                  gp = grid::gpar(cex = size, col = color))
  grid::popViewport()
}


SigLevel<-function(vector){
  return(sapply(vector, function(x) ifelse(x>0.001 && x<0.05, "*", ifelse(x<0.001, "**", "NS"))))
}

FindContrasts<-function(object){
  return(print(names(object[[3]])))
}
