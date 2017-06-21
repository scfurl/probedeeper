


PDVenn<-function(PDobj, index.sel, colors=PDobj@ColObj@match$line, plottype = "ChowRuskey", print_results=F, return_plot=T, return_obj=F){
  #This function takes a PD obj and makes Venn diagrams of the various contrasts
  #index.sel= vector of indices of DE contrasts to plot Venn
  #colors= vector of colors for index.sel
  classvec.sel=PDobj@ColObj@classvec
  VList=list()
  for(i in 1:length(index.sel)){
    VList[[i]]<- rownames(PDobj@LimmaObj@DEGenes[[index.sel[i]]])
  }
  colordriver<-ColObj@match$line[names(ColObj@match$line) %in% PDobj@LimmaObj@Contrasts$comp1[index.sel]]
  names(VList)<-PDobj@LimmaObj@Contrasts$meaning[index.sel]
  Vstem <- Venn(VList)
  Tstem <- compute.Venn(Vstem)
  gp <- VennThemes(Tstem, colourAlgorithm = "sequential", increasingLineWidth=FALSE)
  gps <- gp[["Set"]]
  nSets <- length(gps)
  # substr(names(PD@LimmaObj@DEGenes),1,1)
  # venncolors<-colors[which(LETTERS[seq(from=1, to=length(colors))] %in% substr(names(limma.out[[3]]),1,1)[index.sel])]
  # if(length(venncolors)<nSets){
  #   venncolors<-c(venncolors, RColorBrewer::brewer.pal(12, "Paired")[1:(nSets-length(venncolors))])
  # }
  # make.unique(venncolors)
  for (ix in 1:nSets) {
    gps[[ix]]$col<-colordriver[ix]
  }
  gp[["Set"]] <- gps
  for (ix in 1:nSets) {
    gp$SetText[[ix]]$col<-colordriver[ix]
    gp$SetText[[ix]]$fontsize<-16
  }
  if(print_results){
    for(i in 1:length(VList)){
      print(paste(names(VList[i]), length(VList[[i]]), sep="-"))
    }
  }
  if(return_plot){
    return(plot(Vstem, type=plottype, show = list(SetLabels = TRUE), gp=gp))
  }
  if(return_obj){
    return(Vstem)
  }
}

