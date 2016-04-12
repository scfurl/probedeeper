MAboxplot2<-function(gene, array, classvec, cols=brewer.pal(length(levels(classvec)), "Set3"), reorder=NULL){
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
  stats<-pairwise.t.test(df$y, classvec, p.adjust="fdr")
  array.ind<-as.data.frame(which(stats$p.value < 0.05, arr.ind=T))
  array.ind$row<-array.ind$row+1
  rownames(array.ind)<-NULL
  colnames(array.ind)<-c("start","end")
  array.ind$y<-seq(from=maxh+0.5*(spread), by=spread/2, length.out=nrow(array.ind))
  g<-ggplot(df, aes(x = gp, y = y)) +
    geom_boxplot(size=0.5, alpha=0.6, fill=cols,outlier.size=NULL, width=0.6) +
    geom_point(size=4, colour=cols[classvec], alpha=0.8,position = position_jitter(width = .1)) +
    labs(list(x = "Treatment Group", y = "Expression Value", title=gene)) +
    theme_bw() +
    theme(axis.text=element_text(size=16),
          axis.title.x=element_text(size=20, vjust=0),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.title.y=element_text(size=16, vjust=0.5), plot.title = element_text(vjust = 0, size=20),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  if(nrow(array.ind)==0){print(g)}else{
  for(i in 1:nrow(array.ind)){
    g<-g+geom_segment(aes_string(x = array.ind$start[i], y = array.ind$y[i], xend = array.ind$end[i], yend=array.ind$y[i]), lwd=0.5,arrow = arrow(angle=90, ends="both", length = unit(0.1, "cm")))
  }
  print(g)}
}


strwrap_strip_text = function(p, pad=0.05) {
  # get facet font attributes
  th = theme_get()
  if (length(p$theme) > 0L)
    th = th + p$theme

  require("grid")
  grobs <- ggplotGrob(p)

  # wrap strip x text
  if ((class(p$facet)[1] == "grid" && !is.null(names(p$facet$cols))) ||
        class(p$facet)[1] == "wrap")
  {
    ps = calc_element("strip.text.x", th)[["size"]]
    family = calc_element("strip.text.x", th)[["family"]]
    face = calc_element("strip.text.x", th)[["face"]]

    if (class(p$facet)[1] == "wrap") {
      nm = names(p$facet$facets)
    } else {
      nm = names(p$facet$cols)
    }

    # get number of facet columns
    levs = levels(factor(p$data[[nm]]))
    npanels = length(levs)
    if (class(p$facet)[1] == "wrap") {
      cols = n2mfrow(npanels)[1]
    } else {
      cols = npanels
    }

    # get plot width
    sum = sum(sapply(grobs$width, function(x) convertWidth(x, "in")))
    panels_width = par("din")[1] - sum  # inches
    # determine strwrap width
    panel_width = panels_width / cols
    mx_ind = which.max(nchar(levs))
    char_width = strwidth(levs[mx_ind], units="inches", cex=ps / par("ps"),
                          family=family, font=gpar(fontface=face)$font) /
      nchar(levs[mx_ind])
    width = floor((panel_width - pad)/ char_width)  # characters

    # wrap facet text
    p$data[[nm]] = unlist(lapply(strwrap(p$data[[nm]], width=width,
                                         simplify=FALSE), paste, collapse="\n"))
  }

  if (class(p$facet)[1] == "grid" && !is.null(names(p$facet$rows))) {
    ps = calc_element("strip.text.y", th)[["size"]]
    family = calc_element("strip.text.y", th)[["family"]]
    face = calc_element("strip.text.y", th)[["face"]]

    nm = names(p$facet$rows)

    # get number of facet columns
    levs = levels(factor(p$data[[nm]]))
    rows = length(levs)

    # get plot height
    sum = sum(sapply(grobs$height, function(x) convertWidth(x, "in")))
    panels_height = par("din")[2] - sum  # inches
    # determine strwrap width
    panels_height = panels_height / rows
    mx_ind = which.max(nchar(levs))
    char_height = strwidth(levs[mx_ind], units="inches", cex=ps / par("ps"),
                           family=family, font=gpar(fontface=face)$font) /
      nchar(levs[mx_ind])
    width = floor((panels_height - pad)/ char_height)  # characters

    # wrap facet text
    p$data[[nm]] = unlist(lapply(strwrap(p$data[[nm]], width=width,
                                         simplify=FALSE), paste, collapse="\n"))
  }

  invisible(p)
}

#####3D#######

SFplot3D<-function(A1,A2,A3, classvector, col){
  library(rgl)
pcadata.s<-cbind(A1, A2, A3)
rownames(pcadata.s)<-classvector
colnames(pcadata.s)<-c("PCA1", "PCA2", "PCA3")
pcadata.s<-as.data.frame(pcadata.s)
data.split<-split(pcadata.s, as.factor(rownames(pcadata.s)))
orglspider(as.matrix(data.split[[1]]), as.character(classvector[classvector==levels(classvector)[1]]), col=gcolor[1],size=20.0, label=T)
for(i in 2:length(data.split)){
  orglspider(as.matrix(data.split[[i]]), as.character(classvector[classvector==levels(classvector)[i]]), col=gcolor[i],add=TRUE)
}
for(i in 1:length(data.split)){
  rgl.spheres(as.matrix(data.split[[i]]), col=gcolor[i], radius=6, add=TRUE)
}
#light3d(theta = 0, phi = 15, x = NULL)
rgl.bg( sphere = FALSE, fogtype = "none", color="white",
        back="lines")
grid3d("x", at = NULL, col = "gray", lwd = 1, lty = 1, n = 5)
grid3d("y", at = NULL, col = "gray", lwd = 1, lty = 1, n = 5)
grid3d("z", at = NULL, col = "gray", lwd = 1, lty = 1, n = 5)
}

pdfplot<-function(functionobject, genes, data.sel, classvec, gcolor, filename=paste(prefix, "plot.pdf", sep=""), height=8.5, width=11)
  {
  prefix<-format(Sys.Date(), format="%Y%m%d")
  pdf(filename, height=height, width=width)
  for(i in 1:length(genes)){
    functionobject(genes[i], data.sel, classvec, gcolor)
    print(genes[i])
  }
  dev.off()
}

IsAnnotated<-function(vector, data.mat){
  #This FUnction will return a vector of those probes in a data matrix with rows as genes, samples as columns
  vector<-as.character(vector)
  found<-vector[vector %in% rownames(data.mat)]
  return(found)
}

MAboxplot3<-function(gene, array, limma.obj, classvec, cols=brewer.pal(length(levels(classvec)), "Set3"),
                     reorder=NULL, stat.test="limma", annotate=TRUE){
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
    if(is.null(reorder)==FALSE){stat.df<-stat.df[,reorder]
    }
    rownames(stat.df)<-levels(classvec)[order(levels(classvec))][2:length(levels(classvec))]
    stats<-list(p.value=stat.df)
    array.ind<-as.data.frame(which(stats$p.value < 0.05, arr.ind=T))
    array.ind$row<-match(rownames(stat.df)[array.ind$row], levels(classvec))
    rownames(array.ind)<-NULL
    colnames(array.ind)<-c("start","end")
    array.ind$y<-seq(from=maxh+0.5*(spread), by=spread/2, length.out=nrow(array.ind))
  }
 if(stat.test=="pairwiset"){
   array.ind<-as.data.frame(which(stats$p.value < 0.05, arr.ind=T))
   rownames(array.ind)<-NULL
   array.ind$row<-array.ind$row+1
   colnames(array.ind)<-c("start","end")
   array.ind$y<-seq(from=maxh+0.5*(spread), by=spread/2, length.out=nrow(array.ind))
  }
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


stats.table<-function(df){
  #This function takes a df with three columns, p.value, Comp1 and Comp2 in long format and converts
  # to a df in wide format, eliminating duplicate comparisons
  tmp<-vector()
  tmp2<-vector()
  tmp3<-vector()
  EffComp<-vector()
  for(i in 1:nrow(df)){
    tmp[i]<-c(as.character(df$Comp1[i]))
    tmp2[i]<-c(as.character(df$Comp2[i]))
    tmp3<-c(tmp[i], tmp2[i])
    EffComp[i]<-paste(tmp3[order(tmp3)][1], tmp3[order(tmp3)][2], sep="v")
  }
  df$EffComp<-EffComp
  df$p.value[duplicated(df$EffComp)]<-NA
  df<-unstack(df, p.value~Comp2)
}

# pvalue.lookup<-function(i,j){
#   i<-factor(sapply(limma.master[[1]], function(x) substring(x,1,1)))
#   j<-factor(sapply(limma.master[[1]], function(x) substring(x,2,2)))
#
# }

extractMAT<-function(M, dtype="greater", cutoff=0, replace=0){
  if(dtype=="greater"){
    for(i in 1:nrow(M)){
      for(j in 1:ncol(M)){
        if(M[i,j]<cutoff){
          M[i,j]<-replace
        }
      }
    }
  }
  if(dtype=="less"){
    for(i in 1:nrow(M)){
      for(j in 1:ncol(M)){
        if(M[i,j]>cutoff){
          M[i,j]<-replace
        }
      }
    }
  }
  return(M)
}

quickVenn<-function(list, colors=brewer.pal(length(list),"Paired"), plottype = "ChowRuskey"){
  Vstem <- Venn(list)
  Tstem <- compute.Venn(Vstem)
  gp <- VennThemes(Tstem, colourAlgorithm = "sequential", increasingLineWidth=FALSE)
  gps <- gp[["Set"]]
  nSets <- length(gps)
  venncolors<-colors
  for (ix in 1:nSets) {
    gps[[ix]]$col<-venncolors[ix]
  }
  gp[["Set"]] <- gps
  for (ix in 1:nSets) {
    gp$SetText[[ix]]$col<-venncolors[ix]
    gp$SetText[[ix]]$fontsize<-16
  }
  plot(Vstem, type=plottype, show = list(SetLabels = TRUE), gp=gp)
  print(Vstem)
  return(Vstem)
}
