####GSEA Support####
# MSigDBDir<-"~/Kean Lab/Bioinformatics Resources/MSigDB"
# MSigDBFiles<-list.files(MSigDBDir, full.names=TRUE)
# masterDB<-list()
# for(i in 1:length(MSigDBFiles)){
# mSigDB <- readLines(MSigDBFiles[i])
# mSigDB <- strsplit(mSigDB, "\t")
# names(mSigDB) <- sapply(mSigDB, function(x) x[1])
# mSigDB <- sapply(mSigDB, function(x) x[3:length(x)])
# masterDB<-c(masterDB,mSigDB)}
# save(masterDB, file="~/Kean Lab/Bioinformatics Resources/MSigDB/MasterMsigDB")
# workpath<-"/Bioinformatics Resources"
# path<-paste(rootpath, workpath, sep="")
# setwd(path)
# load("MSigDB/masterDB.150602")
# #library(cnvGSA)

gsea.write.gct <- function(exprMat, gctFn){
  nGenes = nrow(exprMat)
  nConds = ncol(exprMat)
  write("#1.2", file = gctFn, append = F) # All .gct files require this dummy header line.
  write(paste(nGenes, nConds, sep = "\t"), file = gctFn, append = T)
  if(length(make.unique(colnames(exprMat)))!=length(unique(colnames(exprMat))))
  {colnames(exprMat)<-make.unique(colnames(exprMat))}
  write(paste("Name", "Description", paste(colnames(exprMat), collapse = "\t"), sep = "\t"), file = gctFn, append = T)
  # The second column of the .gct file, "Description", is filled out with "na"'s.
  rownames(exprMat) = paste(rownames(exprMat), "na", sep = "\t") # Append "\tna" to every gene name.
  write.table(exprMat, file = gctFn, append = T, quote = F, sep = "\t",
              na = "", row.names = T, col.names = F)
}

gsea.write.cls <- function(classLabels, clsFn, reverse=FALSE){
  nConds = length(classLabels)
  uniqueLabels = levels(as.factor(classLabels))
  if(reverse==TRUE){uniqueLabels <- rev(levels(as.factor(classLabels)))}
  write(paste(nConds, length(uniqueLabels), "1", sep = " "), file = clsFn, append = F)
  write(paste("#", paste(uniqueLabels, collapse=" "), sep=" "), file = clsFn, append = T)
   write(as.character(classLabels), file = clsFn, append = T, sep = " ", ncolumns = nConds)
}

gsea.write.gmt<-function(x, fil){ z <- deparse(substitute(x))
                                  nams=names(x)
                                  cat(nams[1],"NULL" , x[[1]], "\n", sep="\t",
                                      file=fil)
                                  if(length(x)>1){
                                    for (i in 2:length(x) ){ cat(nams[i],"NULL" , x[[i]], "\n",
                                                              file=fil, sep="\t",append=TRUE) }
                                  }
}


readGMTwitheset<-function(file, eset){
  if(file.exists(file) == FALSE) stop("File does not exist")
  tmp<-read.table(file, sep="\t", stringsAsFactors=FALSE)
  tmp<-tmp[3:nrow(tmp),1]
  if(class(eset)=="ExpressionSet"){tmp2<-exprs(eset)}
  else(tmp2<-as.matrix(eset))
  return(tmp2[which(rownames(tmp2) %in% as.character(tmp)),])}

GMTtoHM.list<-function(MSIGName, eset){
  if(class(MSIGName)!="character"){stop("Input is not a character vector")}
  if(class(eset)=="ExpressionSet"){tmp2<-exprs(eset)}
  else(tmp2<-as.matrix(eset))
  output<-list()
  for(i in 1:length(MSIGName)){
    output[[i]]<-rownames(tmp2)[which(rownames(tmp2) %in%
        unlist(masterDB[names(masterDB)==MSIGName[i]],
               use.names=FALSE))]
    names(output)[i]<-MSIGName[i]
  }
  return(output)
}

GMTtoHM.pdf<-function(GMTlist, eset, Fname, cols, sidecols){
  if(class(GMTlist)!="list"){stop("Input is not a list")}
  if(class(eset)=="ExpressionSet"){tmp2<-exprs(eset)}
  else(tmp2<-as.matrix(eset))
  if(class(Fname)!="character"){stop("Filename is not valid")}
  pdf(file= Fname, width=12, height=8)
  for(i in 1:length(GMTlist)){
    HM<-tmp2[which(rownames(tmp2) %in% GMTlist[[i]]),]
    heatmap.2(HM, dendrogram="both", Colv=TRUE, Rowv=TRUE, col=cols,scale="row", trace="none", main=names(GMTlist)[i],
              ColSideColors=sidecols, keysize = 0.9,cexRow=0.6, cexCol=0.6,srtCol=45,density.info="none")

  }
    dev.off()
}

IPAtoGMT<-function(file, output, name){
  if(file.exists(file)==FALSE) {break("File not found")} else{
    tmp<-read.table(file,header = FALSE,fill=TRUE, skip=3, stringsAsFactors=FALSE)
    tmp2<-tmp$V1
    annotation_rhesus<-read.csv("~/Kean Lab/Bioinformatics Data/Rhesus Arrays/Rhesus Annotation/Master.Annotation.csv",colClasses='character',comment.char='#')
    row.names(annotation_rhesus)<-annotation_rhesus[,1]
    annotation_rhesus[,1]<-NULL
    tmp3<-tmp2[match(annotation_rhesus$Gene.Symbol,tmp2)]
    tmp4<-unique(tmp3[!is.na(tmp3)])
    tmp5<-c(name, "NA", tmp4)
    write(tmp5, file = output, ncolumns = length(tmp5), append = FALSE, sep = "\t")
  }}


GCTtoMAT<-function(file){
  if(file.exists(file)==FALSE){stop("File is not valid")}
  MAT<-read.table(file, skip=2, header=TRUE, sep="\t", quote="")
  MAT$Description<-NULL
  names<-MAT$Name
  MAT<-data.matrix(MAT[,2:ncol(MAT)])
  rownames(MAT)<-names
  return(MAT)
}


EnrichmentPlot.SF<-function(filelist, title, plotcols=rep("black", length(filelist)), stats){
  #if(groupnames==NULL){groupnames<-filelist}
  if(class(filelist)!="character"){stop("File List does not appear to be in correct form")}
  if(length(plotcols)!=length(filelist)){stop("Colors not right length")}
  dflist<-list()
  lRES<-vector()
  uRES<-vector()
  for(i in 1:length(filelist)){
    #Rep<-read.delim(filelist[1], header=TRUE, stringsAsFactors=FALSE)
Rep<-read.delim(filelist[i], header=TRUE, stringsAsFactors=FALSE)
df<-cbind.data.frame(as.numeric(Rep$RES), as.numeric(Rep$LIST.LOC), as.character(Rep$CORE_ENRICHMENT))
colnames(df)<-c("RES", "LOC", "ENRICHMENT")
df$max<-df$LOC[which(df$RES == max(df$RES))][1]
df$min<-df$LOC[which(df$RES == min(df$RES))][1]

lRES[i]<-df$RES[which(df$RES == min(df$RES))][1]
uRES[i]<-df$RES[which(df$RES == max(df$RES))][1]
if(abs(lRES[i])>(abs(uRES[i]))){df$dotted<-df$min} else {df$dotted<-df$max}
df$yrug<--(.2+(.1*(i-1)))
dflist[[i]]<-df
      }
ystats<-vector()
for(i in 1:length(filelist)){
    df<-dflist[[i]]
    df$yrug<-min(lRES)-(.2+(.1*(i-1)))
    dflist[[i]]<-df
  }
fdrstats<-as.data.frame(t(round(stats[title,], digits=3)))
  #fdrstats<-as.data.frame(t(round(fdr[title,], digits=3)))
  fdrstats$lRES<-min(lRES)-0.2
colnames(fdrstats)<-c("Name", "lRES")
g<-ggplot() +
  geom_line(data=dflist[[1]], aes(x=LOC, y=RES), colour=plotcols[1])+
  geom_vline(data=dflist[[1]], aes(xintercept=dotted), colour=plotcols[1], linetype="dotted", size=0.8)+
  geom_abline(intercept=0, slope=0)+
  geom_point(data=dflist[[1]],aes(x=LOC, y=yrug, size=as.factor(ENRICHMENT)),colour=plotcols[1], shape=124)+
  #geom_text(data=stats, colour=plotcols[1], aes(label = fdrstats[1], x = -200, y = -0.2))+
  theme_bw()
   for(i in 2:length(dflist)){
  g<-g+geom_line(data=dflist[[i]], aes(x=LOC, y=RES), colour=plotcols[i])
  g<-g+geom_vline(data=dflist[[i]], aes(xintercept=dotted), colour=plotcols[i], linetype="dotted", size=0.8)
  g<-g+geom_point(data=dflist[[i]],aes(x=LOC, y=yrug, size=as.factor(ENRICHMENT)),colour=plotcols[i], shape=124)
  g<-g+ scale_size_manual(values = c(3,6))
  #g<-g+ scale_size_discrete(guide=FALSE)
  g<-g+ theme(legend.position="none")
  g<-g+ labs(title=title)
     }
g <-g + geom_text(data=fdrstats, aes(label=Name, x = rep(-200,length(filelist)), y = seq(min(lRES), (min(lRES)-(length(filelist)-1)*0.1), -0.1)), colour=plotcols, parse=FALSE, size=4)

print(g)
return(g)
}


readGMT.SF<-function(filename){
  if (missing(filename)) {
    stop("Missing 'filename' argument")
  }
  gs.name <- character()
  gs.description <- character()
  gs.genes <- list()
  con <- file(filename, "r")
  recnum <- 0
  repeat {
    hrec <- readLines(con, 1)
    if (length(hrec) == 0) {
      break
    }
    recnum <- recnum + 1
    fields <- unlist(strsplit(hrec, "\t"))
    gs.name[recnum] <- fields[1]
    gs.description[recnum] <- fields[2]
    gs.genes[[recnum]] <- fields[3:length(fields)]
  }
  close(con)
  names(gs.genes)<-gs.name
  return(gs.genes)
}


LEmatrix<-function(LEList, Limma.df, GMTdata, onlyLEgenes=TRUE){
  LEList<-lapply(LEList, function(x) x[order(x)])
  GMTdata<-lapply(LEList, function(x) x[order(x)])
  if(onlyLEgenes==FALSE){pathgenes<-unique(unlist(GMTdata))}
  else{pathgenes<-unique(unlist(LEList))}
  pathgenes<-pathgenes[pathgenes %in% rownames(Limma.df)]
  LFC<-Limma.df[pathgenes,1]
  LFC.df<-cbind(pathgenes,LFC)
  GMT.df<-ldply(GMTdata, rbind)
  LE.df<-ldply(LEList, rbind)
  GMT.df<-GMT.df[order(GMT.df$.id),]
  LE.df<-LE.df[order(LE.df$.id),]
  Output.df<-data.frame()
  for(i in 1:nrow(GMT.df)){
    indices<-which(pathgenes %in% LE.df[i,][!isNA(LE.df[i,])])
    tobind<-as.vector(rep(0, length(pathgenes)))
    for(j in 1:length(indices)){
      tobind[indices[j]]<-as.numeric(LFC.df[indices[j],2])
    }
    Output.df<-rbind(Output.df,tobind)
  }
  rownames(Output.df)<-GMT.df$.id
  colnames(Output.df)<-pathgenes
  return(Output.df)
}

c25 <- c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")

findcommon<-function(pvec, num=1, gmt=masterDB){
  #### THIS FUNCTION WILL TAKE NAMES OF MSIGDB PATHWAYS AND FIND GENES COMMON TO
  ####  THOSE PATHWAYS AND RETURN THE NUMBER OF REPLICATES.  argument pvec = vector of pathway names
  ####  argument num =  the number of pathways to threshold on
  if(class(pvec)!="character"){stop("Data not in correct form")}
    if(num<1){stop("Error in threshold")}
  wlist<-masterDB[which(names(gmt) %in% pvec)]
wtable<-table(as.character(unlist(wlist)))
return(wtable[wtable>=num])
}


###GSEA on Unaffected Genes####

MultipleClassGSEA<-function(data.GSEA, Comparison1, Comparison2, GMTList, classvec=classvec.sel){
datestamp<-format(Sys.Date(), format="%Y%m%d")
rdirectory<-paste("GSEA/", datestamp, "/", sep="")
dir.create(rdirectory, showWarnings=FALSE)
fdirectory<-paste(rdirectory, "Files-", datestamp, "/", sep="")
dir.create(fdirectory, showWarnings=FALSE)
odirectory<-paste(rdirectory, "Output-", datestamp, "/", sep="")
dir.create(odirectory, showWarnings = FALSE)
GMTfile<-paste(fdirectory, "GMTFile.gmt", sep="")
gsea.write.gmt(GMTList, GMTfile)
EnrichmentData<-list()
GSEAcomplist<-list(Number=1:length(Comparison1), Comp1=Comparison1, Comp2=rep(Comparison2, length(Comparison1)),
                   Prefix=paste(str_replace_all(str_replace_all(Comparison1, c("[[:lower:]]"), ""), "([.])", ""),
                                str_replace_all(str_replace_all(Comparison2, c("[[:lower:]]"), ""), "([.])", ""), sep="v"),
                   Signreverse=!substring(Comparison1,1,1)<substring(Comparison2,1,1))
for(i in 1:length(GSEAcomplist[["Number"]])){
  file.cls<-paste(fdirectory, GSEAcomplist$Prefix[i], ".cls", sep="")
  data.cls<-as.character(classvec[classvec %in% c(GSEAcomplist$Comp1[i],GSEAcomplist$Comp2[i])])
  gsea.write.cls(data.cls, file.cls)
  file.gct<-paste(fdirectory, GSEAcomplist$Prefix[i], ".gct", sep="")
  col.sel<-classvec %in% c(GSEAcomplist$Comp1[i],GSEAcomplist$Comp2[i])
  data.gct<-data.GSEA[,col.sel]
  gsea.write.gct(data.gct, file.gct)
  wdirectory<-paste(odirectory, GSEAcomplist$Prefix[i], "/",sep="")
  dir.create(wdirectory, showWarnings = FALSE)
  RunGSEA(paste(fdirectory, GSEAcomplist$Prefix[i], ".gct", sep=""),
          paste(fdirectory, GSEAcomplist$Prefix[i], ".cls", sep=""),
          GMTfile,
          wdirectory,
          GSEAcomplist$Prefix[i],
          GSEAcomplist$Signreverse[i])
  if(EnrichmentPlot==TRUE){
    filelist<-list()
    for(j in 1:length(GMTList)){
      regex<-paste(GSEAcomplist$Prefix[i], names(GMTList)[j], sep=".")
      tmp2<-grep(regex, list.files(wdirectory), value=TRUE)
      filelist[[i]][j]<-0
        #EnrichmentData[[j]][[i]]<-NULL
    }
  }
}
}


gsea.read.gct<-function(fn="NULL", allowDuplicatedRows=FALSE){
  df<-read.delim(fn, skip=2, stringsAsFactors=FALSE)
  df<-df[-which(colnames(df)=="Description")]
  if(allowDuplicatedRows==TRUE){
    rownames(df)<-make.unique(df$Name)
    df<-df[-which(colnames(df)=="Name")]
    return(as.matrix(df))
  }
  else
  {
    df<-df[!duplicated(df$Name),]
    rownames(df)<-df$Name
    warning("Duplicates Removed")
    df<-df[-which(colnames(df)=="Name")]
    return(as.matrix(df))
  }
}



readGCTasDT<-function(fn="NULL", allowDuplicatedRows=FALSE, designator="Name", description_col_blank=F, description_designator="Description"){
  if(description_col_blank){
    skip<-3
    cn<-unlist(strsplit(readLines(file.path(RES_DIR, "ssGSEA", "180713", "-combined.gct"), n=3)[3], "\t"))
    colnames<-cn[-c(1:2)]
    designator<-cn[1]
    df<-data.table::fread(fn, skip=skip, stringsAsFactors=FALSE)
    colnames(df)<-c(designator, colnames)
    df[[description_designator]]<-rep("na", nrow(df))
    data.table::setcolorder(df, c(designator, description_designator, colnames[-1]))
  }else{
    skip<-2
    df<-data.table::fread(fn, skip=skip, stringsAsFactors=FALSE)
    df<-df[-which(colnames(df)==description_designator)]
  }
  if(allowDuplicatedRows==TRUE){
    rownames(df)<-make.unique(df[[designator]])
    df<-df[-which(colnames(df)==designator)]
    return(df)
  }
  else{
    df<-df[!duplicated(df[[designator]]),]
    rownames(df)<-df[[designator]]
    warning("Duplicates Removed")
    df<-df[-which(colnames(df)==designator)]

    return(df)
  }
}

gsea.heatmap.matrix<-function(V){
    # Plots a heatmap "pinkogram" of a gene expression matrix

    n.rows <- length(V[,1])
    n.cols <- length(V[1,])
    row.mean <- apply(V, MARGIN=1, FUN=mean)
    row.sd <- apply(V, MARGIN=1, FUN=sd)
    row.n <- length(V[,1])
    for (i in 1:n.rows) {
      if (row.sd[i] == 0) {
        V[i,] <- 0
      } else {
        V[i,] <- (V[i,] - row.mean[i])/(0.5 * row.sd[i])
      }
      V[i,] <- ifelse(V[i,] < -6, -6, V[i,])
      V[i,] <- ifelse(V[i,] > 6, 6, V[i,])
    }
    heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
    heatm[1:n.rows,] <- V[seq(1, n.rows, 1),]
    rownames(heatm)<-rownames(V)
    colnames(heatm)<-colnames(V)
    return(heatm)
}

clusterToGMT<-function(kmeans, moduleNames=NULL){
  if(class(kmeans)!="kmeans"){
    stop("input not of correct class")
  }
  else{
    output<-list()
    for(i in range(kmeans$cluster)[1]:range(kmeans$cluster)[2]){
      output[[i]]<-names(kmeans$cluster[kmeans$cluster==i])
    }
    if(length(moduleNames)==0){
      names(output)<-levels(as.factor(kmeans$cluster))
    }
    else
    {
      names(output)<-moduleNames
    }
  }
  return(output)
}

GMTtoList<-function(GMTfn){
  file.data<-readLines(GMTfn)
  num.lines<-length(file.data)
  output<-vector("list", num.lines)
  for(i in 1:num.lines){
    vec<-unlist(strsplit(file.data[i], "\t"))
    output[[i]]<-vec[3:length(vec)]
    names(output)[i]<-vec[1]
  }
  return(output)
}

read.odfdata<-function(odf.filename){
  file.data<-readLines(odf.filename)
  num.lines<-length(file.data)
  headerlines<-as.numeric(strsplit(file.data[2], " ")[[1]][2])
  header.list<-vector("list", headerlines+2)
  for(i in 1:(headerlines+2)){
    vec<-unlist(strsplit(file.data[i], "\t"))
    if(length(vec)==1){
      vec2<-unlist(strsplit(vec, " "))
      header.list[[i]]<-vec2[2]
      names(header.list)[i]<-vec2[1]
      }
    else{
        header.list[[i]]<-vec[2:length(vec)]
        names(header.list)[i]<-vec[1]
      }
    }
  df<-NULL
  for(i in (headerlines+3):num.lines){
    df<-rbind(df,unlist(strsplit(file.data[i], "\t")))
  }
  type.ind<-which(names(header.list)=="COLUMN_TYPES:")
  name.ind<-which(names(header.list)=="COLUMN_NAMES:")
  df2<-assignTypeDataFrame(as.data.frame(df, stringsAsFactors=FALSE), header.list[[type.ind]])
  colnames(df2)<-header.list[[name.ind]]
  return(df2)
}

assignTypeDataFrame<-function(df, coltypes){
  coltypes.r<-matchTypes(coltypes)
  num.int<-which(coltypes.r %in% "numeric")
  char.int<-which(coltypes.r %in% "character")
  df[,num.int]<-sapply(df[,num.int], as.numeric)
  df[,char.int]<-sapply(df[,char.int], as.character)
  return(df)
}

matchTypes<-function(vector){
  vector.fac<-as.factor(tolower(vector))
  check<-levels(vector.fac)
  allowed<-c("character", "numeric", "float", "int", "string")
  if(length(which(!check %in% allowed)==0)){break("input type not allowed")}
  vector<-tolower(vector)
  vector[vector=="float"]<-"numeric"
  vector[vector=="int"]<-"numeric"
  vector[vector=="string"]<-"character"
  return(vector)
}

setClass("GenePatternObject", representation(
  Contents="list",
  #Contains content data
  ODFHeader="list",
  #Contains metadata
  ODFData="data.frame",
  GCTData="data.frame"))

newGenePatternObject<-function(odf.fn, gct.fn){
  odf<-read.odf(odf.fn)
  Contents<-list(ODFFileName=odf.fn, GCTFileName=gct.fn)
  newGPO<-new("GenePatternObject", Contents=Contents, ODFHeader=odf[[1]],
              ODFData=odf[[2]], GCTData=as.data.frame(gsea.read.gct(gct.fn, allowDuplicatedRows=TRUE)))
}

read.odf<-function(odf.filename){
  file.data<-readLines(odf.filename)
  num.lines<-length(file.data)
  headerlines<-as.numeric(strsplit(file.data[2], " ")[[1]][2])
  header.list<-vector("list", headerlines+2)
  output<-list()
  for(i in 1:(headerlines+2)){
    vec<-unlist(strsplit(file.data[i], "\t"))
    if(length(vec)==1){
      vec2<-unlist(strsplit(vec, " "))
      header.list[[i]]<-vec2[2]
      names(header.list)[i]<-vec2[1]
    }
    else{
      header.list[[i]]<-vec[2:length(vec)]
      names(header.list)[i]<-vec[1]
    }
  }
  output[[1]]<-header.list
  df<-NULL
  for(i in (headerlines+3):num.lines){
    df<-rbind(df,unlist(strsplit(file.data[i], "\t")))
  }
  type.ind<-which(names(header.list)=="COLUMN_TYPES:")
  name.ind<-which(names(header.list)=="COLUMN_NAMES:")
  df2<-assignTypeDataFrame(as.data.frame(df, stringsAsFactors=FALSE), header.list[[type.ind]])
  colnames(df2)<-header.list[[name.ind]]
  output[[2]]<-df2
  names(output)<-c("ODFHeader", "ODFData")
  return(output)
}
