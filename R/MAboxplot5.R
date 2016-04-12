#####MABoxPlot6#####
#MABoxPlot4 
#Bug fix to not plot NA statistical values as significant
#New feature to allow input of statistical cutoff
#150809 Fixed bug that prevented boxplot of only two groups.
#MABoxPlot5
#New feature to plot sample names on plot
#MABoxPlot6
#Allow for further manipulation of color

# #Dummy Data
# gene<-AURKA.PROBES[2]
# sampleNames<-paste(finaleset$sampleNames, finaleset$Batch, sep="-")
# array<-data
# limma.obj<-limma.out
# classvec<-classvec
# stat.test="pairwiset"
# p.value=0.05

MAboxplot5<-function(gene, array, limma.obj, classvec, cols=brewer.pal(length(levels(classvec)), "Set3"), 
                     reorder=NULL, stat.test="limma", annotate=TRUE, p.value=0.05, sampleNames=NULL){
  library(grid)
  classvec<-as.factor(classvec)
  if(is.null(reorder)==FALSE){classvec<-factor(classvec,levels(classvec)[reorder])
                              cols<-cols[reorder]  
  }
  obj<-array[gene,]
  df<-data.frame(gp=classvec, y=obj)
  maxh <- max(df$y)
  minh<-min(df$y)
  spread<-(maxh-minh)/14
  maxh<-maxh+spread
  df$maxh<-maxh
  df$spread<-spread
  ds <- ddply(df, .(gp), summarise, mean = mean(y), sd = sd(y))
  if(stat.test=="pairwiset"){
    stats<-pairwise.t.test(df$y, classvec, p.adjust="fdr")
  }
  if(stat.test=="limma"){
    stat.df<-data.frame(Comp1=as.factor(limma.obj[[5]]$Comp1), Comp2=as.factor(limma.obj[[5]]$Comp2), p.value=ExtractLIMMA(limma.obj, gene)$adj.P.Val, stringsAsFactors=FALSE)
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
    g<-ggplot(df, aes(x = gp, y = y)) +
      geom_boxplot(size=0.5, alpha=0.6, fill=cols,outlier.size=NULL, width=0.6) +
      geom_point(size=4, colour=cols[classvec], alpha=0.8,position = position_jitter(width = .1)) +
      geom_text(data=df, aes(x = gp, y = y, label=sampleNames), size = 3, hjust=-1) +
      labs(list(x = "Treatment Group", y = "Log2 Transformed Data", title=gene)) +
      #ylab(expression(paste("Log", [2], " Transformed Data", sep="")))+
      theme_bw() +
      theme(axis.text=element_text(size=16),
          axis.title.x=element_text(size=20, vjust=1),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.title.y=element_text(size=16, vjust=0.5), plot.title = element_text(vjust = 0, size=20),
          axis.line = element_line(colour = "black"),
          #text=element_text(family="Myriad Pro"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
    if(nrow(array.ind)==0){
      print(g)
    }
    else
    {
      for(i in 1:nrow(array.ind)){
        g<-g+geom_segment(aes_string(x = array.ind$start[i], y = array.ind$y[i], xend = array.ind$end[i], yend=array.ind$y[i]), lwd=0.5,arrow = arrow(angle=90, ends="both", length = unit(0.1, "cm")))
      }
    }
    print(g)
    footie1<-ifelse(stat.test=="limma", paste("Mod. Bayesian T statistic corrected using ", limma.obj[[6]]$p.adjust, sep=""), "Pairwise T Test, FDR-corrected")
    Footnote.txt<-paste("Horizontal bars indicate p <0.05 using ", footie1, sep="")
    makeFootnote(Footnote.txt,  color = "black")
  }
  else
  {
      g<-ggplot(df, aes(x = gp, y = y)) +
        geom_boxplot(size=0.5, alpha=0.6, fill=cols,outlier.size=NULL, width=0.6) +
        geom_point(size=4, colour=cols[classvec], alpha=0.8,position = position_jitter(width = .1)) +
        labs(list(x = "Treatment Group", y = "Log2 Transformed Data", title=gene)) +
        #ylab(expression(paste("Log", [2], " Transformed Data", sep="")))+
        theme_bw() +
        theme(axis.text=element_text(size=16),
              axis.title.x=element_text(size=20, vjust=1),
              axis.text.x = element_text(angle = 45, hjust=1),
              axis.title.y=element_text(size=16, vjust=0.5), plot.title = element_text(vjust = 0, size=20),
              axis.line = element_line(colour = "black"),
              #text=element_text(family="Myriad Pro"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())
      if(nrow(array.ind)==0){
        print(g)
      }
      else
      {
        for(i in 1:nrow(array.ind)){
          g<-g+geom_segment(aes_string(x = array.ind$start[i], y = array.ind$y[i], xend = array.ind$end[i], yend=array.ind$y[i]), lwd=0.5,arrow = arrow(angle=90, ends="both", length = unit(0.1, "cm")))
        }
      }
      print(g)
      footie1<-ifelse(stat.test=="limma", paste("Mod. Bayesian T statistic corrected using ", limma.obj[[6]]$p.adjust, sep=""), "Pairwise T Test, FDR-corrected")
      Footnote.txt<-paste("Horizontal bars indicate p <0.05 using ", footie1, sep="")
      makeFootnote(Footnote.txt,  color = "black")
    }
}

# write a simple function to add footnote
makeFootnote <- function(footnoteText =
                           format(Sys.time(), "%d %b %Y"),
                         size = .7, color = grey(.5))
{
  require(grid)
  pushViewport(viewport())
  grid.text(label = footnoteText ,
            x = unit(1,"npc") - unit(2, "mm"),
            y = unit(2, "mm"),
            just = c("right", "bottom"),
            gp = gpar(cex = size, col = color))
  popViewport()
}