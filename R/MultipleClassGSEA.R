

setClass("MultipleClassGSEAObject", representation(
  ObjectInfo="list",
    #Contains GSEAcomplist, GMTList, timestamp, odirectory
  Data="list",
    #Contains df of Comparisons grouped by GMTList
  PlotData="list", Stats="list"))
    #Contains statistical data for Enrichment plots


###GSEA on Unaffected Genes####
MultipleClassGSEA<-function(data.GSEA, Comparison1, Comparison2, GMTList, classvec=NULL, runGSEAcode=TRUE, reshuffling.type, directory=getwd(), uniquelabel=""){
#   os<-Sys.info()['sysname']
#   if(os=="Darwin"){GSEA.program.location <- "/Library/Frameworks/R.framework/Resources/library/GSEA/GSEA.1.0.R"}   #  R source program (change pathname to the rigth location in local machine)
#   source(GSEA.program.location, verbose=T, max.deparse.length=9999)
  datestamp<-format(Sys.Date(), format="%Y%m%d")
  if(uniquelabel!=""){
    uniquelabel<-paste("-", uniquelabel, sep="")
  }
  mdirectory<-paste(file.path(directory, "GSEA"))
   rdirectory<-file.path(mdirectory, paste(datestamp, "-", uniquelabel, sep=""))
   fdirectory<-file.path(rdirectory, paste("Files-", datestamp, sep=""))
   odirectory<-file.path(rdirectory, paste("Output-", datestamp, sep=""))

   GMTfile<-file.path(fdirectory, "GMTFile.gmt")
   Comp1<-Comparison1
   Comp2<-Comparison2
   if(length(Comp1)==1){Comp1<-rep(Comparison1, length(Comparison2))}
   if(length(Comp2)==1){Comp2<-rep(Comparison2, length(Comparison1))}
   if(length(Comparison2)!=1 && length(Comparison1)!=1){stop("Input Problem")}
   GSEAcomplist<-list(Number=max(c(length(Comparison1), length(Comparison2))), Comp1=Comp1, Comp2=Comp2,
                      Prefix=paste(stringr::str_replace_all(stringr::str_replace_all(Comparison1, c("[[:lower:]]"), ""), "([.])", ""),
                                   stringr::str_replace_all(stringr::str_replace_all(Comparison2, c("[[:lower:]]"), ""), "([.])", ""), sep="v"),
                      Signreverse=!substring(Comparison1,1,1)<substring(Comparison2,1,1))
   Objectinfo<-list(GSEACompList=GSEAcomplist, GMTList=GMTList, Directory=odirectory, TimeStamp=Sys.time())
  MCGO<-new("MultipleClassGSEAObject", ObjectInfo=Objectinfo,
                            Data=list(length=GSEAcomplist[["Number"]]),
                            PlotData=list(length=GSEAcomplist[["Number"]]))
  MCGO@Data<-rep(list(list()), length(GMTList))
  MCGO@PlotData<-rep(list(data.frame(max=numeric(0), min=numeric(0), lRES=numeric(0), uRES=numeric(0), dotted=numeric(0))), length(GMTList))
  names(MCGO@PlotData)<-names(GMTList)
  statsOut.df<-data.frame()
  stats.df<-list()
  if(runGSEAcode==TRUE){
#     os<-Sys.info()['sysname']
#     if(os=="Darwin"){GSEA.program.location <- "/Library/Frameworks/R.framework/Resources/library/GSEA/GSEA.1.0.R"}   #  R source program (change pathname to the rigth location in local machine)
#     source(GSEA.program.location, verbose=T, max.deparse.length=9999)
    dir.create(rdirectory, showWarnings=FALSE)
    dir.create(fdirectory, showWarnings=FALSE)
    dir.create(odirectory, showWarnings=FALSE)
    gsea.write.gmt(GMTList, GMTfile)
  }
  for(i in 1:GSEAcomplist[["Number"]]){
    file.cls<-file.path(fdirectory, paste(GSEAcomplist$Prefix[i], ".cls", sep=""))
    data.cls<-as.character(classvec[classvec %in% c(GSEAcomplist$Comp1[i],GSEAcomplist$Comp2[i])])
    file.gct<-file.path(fdirectory, paste(GSEAcomplist$Prefix[i], ".gct", sep=""))
    col.sel<-classvec %in% c(GSEAcomplist$Comp1[i],GSEAcomplist$Comp2[i])
    data.gct<-data.GSEA[,col.sel]
    wdirectory<-paste(odirectory, GSEAcomplist$Prefix[i], "/",sep="")
    if(runGSEAcode==TRUE){
    gsea.write.cls(data.cls, file.cls)
    gsea.write.gct(data.gct, file.gct)
    dir.create(wdirectory, showWarnings = FALSE)
    RunGSEA(in.data = file.path(fdirectory, paste(GSEAcomplist$Prefix[i], ".gct", sep="")),
            in.cls =file.path(fdirectory, paste(GSEAcomplist$Prefix[i], ".cls", sep="")),
              in.gmt =  GMTfile,
              directory      = wdirectory,
              prefix            = GSEAcomplist$Prefix[i],
              reshuffling.type      = reshuffling.type,
              signreverse          = GSEAcomplist$Signreverse[i])
    }
      for(j in 1:length(GMTList)){
        ###CREATE DF of MCGO DATA
        regex<-paste("^",GSEAcomplist$Prefix[i], ".", names(GMTList)[j], ".report.",sep="")
        filename<-grep(regex, list.files(wdirectory), value=TRUE)
        Rep<-read.delim(paste(wdirectory, filename, sep=""), header=TRUE, stringsAsFactors=FALSE)
        df<-data.frame()
        #yrug<--(.2+(.1*(i-1)))
        df<-cbind.data.frame(as.numeric(Rep$RES), as.numeric(Rep$LIST.LOC), as.character(Rep$CORE_ENRICHMENT), as.character(Rep$GENE))
        colnames(df)<-c("RES", "LOC", "ENRICHMENT", "GENE")
        MCGO@Data[[j]][[i]]<-df
        names(MCGO@Data[[j]])[i]<-GSEAcomplist$Prefix[i]
        ###CREATE DF of MCGO PLOTDATA
        max=df$LOC[which(df$RES == max(df$RES))][1]
        min=df$LOC[which(df$RES == min(df$RES))][1]
        lRES<-df$RES[which(df$RES == min(df$RES))][1]
        uRES<-df$RES[which(df$RES == max(df$RES))][1]
        if(abs(lRES)>(abs(uRES))){dotted<-min} else {dotted<-max}
        MCGO@PlotData[[j]][i,]<-c(max,min,lRES,uRES, dotted)
        #ystats<-vector()
#         for(i in 1:length(filelist)){
#           df<-dflist[[i]]
#           df$yrug<-min(lRES)-(.2+(.1*(i-1)))
#           dflist[[i]]<-df
#         }
    }
    names(MCGO@Data)<-names(GMTList)
    regex<-paste("^",GSEAcomplist$Prefix[i], ".SUMMARY.RESULTS.REPORT",sep="")
    filename<-grep(regex, list.files(wdirectory), value=TRUE)
    if(length(filename)>1){
      tmp<-read.delim(paste(wdirectory, filename[1], sep=""), header=TRUE, stringsAsFactors=FALSE)
      for(j in 2:length(filename)){
        tmp2<-read.delim(paste(wdirectory, filename[j], sep=""), header=TRUE, stringsAsFactors=FALSE)
        tmp<-rbind(tmp, tmp2[1:nrow(tmp2),])
      }
      stats.df[[i]]<-tmp
    }
    else
    {
    stats.df[[i]]<-read.delim(paste(wdirectory, filename, sep=""), header=TRUE, stringsAsFactors=FALSE)
    }
#     for(j in 1:length(GMTList)){
#       statsOut.df[i,j]<-stats.df$FDR.q.val[j]
#     }
#     for(j in 1:length(GMTList)){
#       statsOut.df[i,j]<-stats.df[[i]]$FDR.q.val[j]
#     }
  }
  MCGO@Stats<-stats.df
  names(MCGO@Stats)<-GSEAcomplist$Prefix
  print(MCGO@ObjectInfo$GSEACompList)
  print(names(MCGO@ObjectInfo[[2]]))
  print(MCGO@ObjectInfo$Directory)
  print(MCGO@ObjectInfo$TimeStamp)
return(MCGO)
}

PlotMultipleClassGSEAObject<-function(MCGO, index, plotcols=rep("Black", length(MCGO@Data[[index]]))){
  dflist<-MCGO@Data[[index]]
  plotdf<-MCGO@PlotData[[index]]
  #statindex<-match(names(MCGO@ObjectInfo[[2]])[index], names(MCGO@ObjectInfo[[2]])[order(names(MCGO@ObjectInfo[[2]]))])
  fdr<-vector()
  for(i in 1:length(MCGO@Data[[index]])){
      statindex<-match(names(MCGO@ObjectInfo[[2]])[index], MCGO@Stats[[i]]$GS)
      fdr[i]<-round(MCGO@Stats[[i]]$FDR.q.val[statindex], 3)
  }
  #fdrstats<-data.frame(FDR=fdr, lRES=rep(min(MCGO@PlotData[[index]]$lRES)-0.2, length(MCGO@Data[[index]])))
  x<-vector()
  y<-vector()
  lRES<-vector()
  lRES=rep(min(MCGO@PlotData[[index]]$lRES), length(MCGO@Data[[index]]))
  for(i in 1:length(fdr)){
    x[i]<--300
    y[i]<-lRES[i]-((i-1)*.1)
  }
  for(i in 1:length(MCGO@Data[[index]])){
    YRUG<-rep(y[i], nrow(dflist[[i]]))
    dflist[[i]]<-cbind(dflist[[i]], YRUG, stringsAsFactors=FALSE)
  }
  fdrstats<-data.frame(FDR=fdr, lRES=lRES, x=x, y=y)
g<-ggplot() +
  geom_line(data=dflist[[1]], aes(x=LOC, y=RES), colour=plotcols[1])+
  geom_vline(data=plotdf[1,], aes(xintercept=dotted), colour=plotcols[1], linetype="dotted", size=0.8)+
  geom_abline(intercept=0, slope=0)+
  geom_point(data=dflist[[1]],aes(x=LOC, y=YRUG, size=as.factor(ENRICHMENT)),colour=plotcols[1], shape=124)+
  #geom_text(data=stats, colour=plotcols[1], aes(label = fdrstats[1], x = -200, y = -0.2))+
  theme_bw()
if(length(dflist)>1){
for(i in 2:length(dflist)){
  g<-g+geom_line(data=dflist[[i]], aes(x=LOC, y=RES), colour=plotcols[i])
  g<-g+geom_vline(data=plotdf[i,], aes(xintercept=dotted), colour=plotcols[i], linetype="dotted", size=0.8)
  g<-g+geom_point(data=dflist[[i]],aes(x=LOC, y=YRUG, size=as.factor(ENRICHMENT)),colour=plotcols[i], shape=124)
  g<-g+ scale_size_manual(values = c(3,6))
  #g<-g+ scale_size_discrete(guide=FALSE)
  g<-g+ theme(legend.position="none")
}
}
g<-g+ labs(title=names(MCGO@Data)[index])
#g <-g + geom_text(data=fdrstats, aes(label=FDR, x = rep(-200, nrow(fdrstats)), y = seq(min(lRES), (min(lRES)-(nrow(fdrstats)-1)*0.1), -0.1)), colour=plotcols, parse=FALSE, size=4)
g <-g + geom_text(data=fdrstats, aes(label=FDR, x=x, y=y), colour=plotcols, parse=FALSE, size=4)
print(g)
return(g)
}

MCGOpdf<-function(MCGO, filename="", plotcols=rep("black", length(names(MCGO@ObjectInfo$GMTList))), gcolor=gcolor){
  pdf(filename, height=8.5, width=11)
  pie(rep(1,length(gcolor)), col=gcolor, labels=names(gcolor))
  for(i in 1:length(names(MCGO@ObjectInfo$GMTList))){
    PlotMultipleClassGSEAObject(MCGO,i, plotcols=plotcols)
  }
  dev.off()
}

RunGSEA<-function(in.data, in.cls, in.gmt, directory, prefix, signreverse, reshuffling.type){
  # GSEA 1.0 -- Gene Set Enrichment Analysis / Broad Institute
  # R script to run GSEA
  ######command###########
  GSEA( # Input/Output Files :-------------------------------------------
        input.ds =  in.data,                 # Input gene expression Affy dataset file in RES or GCT format
        input.cls = in.cls,                 # Input class vector (phenotype) file in CLS format
        gs.db =     in.gmt,          # Gene set database in GMT format
        output.directory      = directory,              # Directory where to store output and results (default: "")
        #  Program parameters :-------------------------------------------------------------------------------------------------------------------------
        doc.string            = prefix,        # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
        non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
        reshuffling.type      = reshuffling.type, # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
        nperm                 = 1000,            # Number of random permutations (default: 1000)
        weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
        nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
        fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
        fdr.q.val.threshold   = 1,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
        topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
        adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
        gs.size.threshold.min = 5,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
        gs.size.threshold.max = 1500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
        reverse.sign          = signreverse,      # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
        preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
        random.seed           = 760435,          # Random number generator seed. (default: 123456)
        perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
        fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
        replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
        save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
        OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
        use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
  )

}

ExtractMCGO<-function(MCGO, name){
  if(length(which(names(MCGO@ObjectInfo$GMTList) %in% name))==0){
   stop("Data Not Found")
  }
  else
  index<-which(names(MCGO@ObjectInfo$GMTList) %in% name)
  dflist<-MCGO@Data[[index]]
  return(dflist)
}

QuickGSEA<-function(data.GSEA, Comp1, Comp2, classvec, GMT, directory, prefix, signreverse=FALSE, reshuffling.type="gene.labels"){
file.cls<-paste(directory, "/", format(Sys.time(), "%a-%b-%d-%Y-%X"),  ".cls", sep="")
data.cls<-as.character(classvec[classvec %in% c(Comp1, Comp2)])
file.gct<-paste(directory, "/", format(Sys.time(), "%a-%b-%d-%Y-%X"),  ".gct", sep="")
col.sel<-classvec %in% c(Comp1, Comp2)
GMTfile<-paste(directory, "/", format(Sys.time(), "%a-%b-%d-%Y-%X"),  ".gmt", sep="")
data.gct<-data.GSEA[,col.sel]
odirectory<-paste(directory, "/Output", format(Sys.time(), "%a-%b-%d-%Y-%X"), "/", sep="")
dir.create(odirectory, showWarnings=FALSE)
  gsea.write.gmt(GMT, GMTfile)
  gsea.write.cls(data.cls, file.cls)
  gsea.write.gct(data.gct, file.gct)
RunGSEA(file.gct, file.cls, GMTfile, odirectory, prefix=format(Sys.time(), "%a-%b-%d-%Y-%X"), signreverse, reshuffling.type)
}

#debug(PlotMultipleEnrichmentPlots)
PlotMultipleEnrichmentPlots<-function(filelist, plotcols=rep("Black", length(filelist)), stats=FALSE){
  dflist<-list()
  PlotData<-data.frame(max=numeric(0), min=numeric(0), lRES=numeric(0), uRES=numeric(0), dotted=numeric(0))
  for(i in 1:length(filelist)){
    dflist[[i]]<-read.delim(filelist[[i]], header=TRUE)
    max=dflist[[i]]$LIST.LOC[which(dflist[[i]]$RES == max(dflist[[i]]$RES))]
    min=dflist[[i]]$LIST.LOC[which(dflist[[i]]$RES == min(dflist[[i]]$RES))]
    lRES<-dflist[[i]]$RES[which(dflist[[i]]$RES == min(dflist[[i]]$RES))]
    uRES<-dflist[[i]]$RES[which(dflist[[i]]$RES == max(dflist[[i]]$RES))]
    if(abs(lRES)>(abs(uRES))){dotted<-min} else {dotted<-max}
    PlotData[i,]<-rbind(c(max,min,lRES,uRES, dotted), PlotData)
    }
  x<-vector()
  y<-vector()
  lRES.global<-vector()
  lRES.global=rep(min(PlotData$lRES), length(filelist))
  for(i in 1:length(filelist)){
    x[i]<--300
    y[i]<-lRES.global[i]-((i-1)*.1)
  }
  for(i in 1:length(filelist)){
    YRUG<-rep(y[i], nrow(dflist[[i]]))
    dflist[[i]]<-cbind(dflist[[i]], YRUG, stringsAsFactors=FALSE)
  }

  g<-ggplot() +
    geom_line(data=dflist[[1]], aes(x=LIST.LOC, y=RES), colour=plotcols[1])+
    geom_vline(data=PlotData[1,], aes(xintercept=dotted), colour=plotcols[1], linetype="dotted", size=0.8)+
    geom_abline(intercept=0, slope=0)+
    geom_point(data=dflist[[1]],aes(x=LIST.LOC, y=YRUG, size=as.factor(CORE_ENRICHMENT)),colour=plotcols[1], shape=124)+
    #geom_text(data=stats, colour=plotcols[1], aes(label = fdrstats[1], x = -200, y = -0.2))+
    theme_bw()
  if(length(filelist)>1){
    for(i in 2:length(filelist)){
      g<-g+geom_line(data=dflist[[i]], aes(x=LIST.LOC, y=RES), colour=plotcols[i])
      g<-g+geom_vline(data=PlotData[i,], aes(xintercept=dotted), colour=plotcols[i], linetype="dotted", size=0.8)
      g<-g+geom_point(data=dflist[[i]],aes(x=LIST.LOC, y=YRUG, size=as.factor(CORE_ENRICHMENT)),colour=plotcols[i], shape=124)
      g<-g+ scale_size_manual(values = c(3,6))
      #g<-g+ scale_size_discrete(guide=FALSE)
      g<-g+ theme(legend.position="none")
    }
  }
  #print(g)
  return(g)
}


ExtractStats<-function(MCGO, index){
  return(MCGO@Stats[[index]][,c(1,2,5,6,7)])
}
