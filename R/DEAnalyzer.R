####DEAnalyzer####

setClass("DEAnalyzerInput", representation(
  InputGeneSet="list",
  DEData="list"))

setClass("DEAnalyzerCalcs", representation(
  ObjectInfo="list",
  InputGeneSet="list",
  DEData="list",
  FoldChange="list",
  Pvalues="list",
  UpDnClassvec="list",
  Up="list",
  Down="list",
  GO.BP.UP="list",
  GO.BP.DN="list",
  GO.MF.UP="list",
  GO.MF.DN="list",
  PW.UP="list",
  PW.DN="list"))


autoLIMMA2DEAnalyzer<-function(limmalist, index.sel){
  InputGeneLists<-list()
  for(i in 1:length(index.sel)){
    InputGeneLists[[i]]<-rownames(limmalist[[3]][[index.sel[i]]])
    names(InputGeneLists)[i]<-limmalist[[2]][[index.sel[i]]][1]
  }
  new("DEAnalyzerInput", InputGeneSet=InputGeneLists,
      DEData=limmalist[[4]][index.sel])
}


Input2DEAnalyzer<-function(InputGeneLists, DEData){
  new("DEAnalyzerInput", InputGeneSet=InputGeneLists,
                                   DEData=DEData)
}

CalcDEAnalyzer<-function(object, GOopt=TRUE, PWopt=TRUE){
  if(class(object)!="DEAnalyzerInput")stop("Input problem")
  calc.object<-new("DEAnalyzerCalcs")
  calc.object@InputGeneSet<-object@InputGeneSet
  calc.object@DEData<-object@DEData
  calc.object@ObjectInfo<-list(NumberGeneSets<-length(object@InputGeneSet),
                                ContainsGOData=GOopt, ContainsPWData=PWopt)
  for(i in 1:length(object@InputGeneSet)){
    FC<-object@DEData[[i]]$logFC[rownames(object@DEData[[i]]) %in% object@InputGeneSet[[i]]]
    PVal<-object@DEData[[i]]$adj.P.Val[rownames(object@DEData[[i]]) %in% object@InputGeneSet[[i]]]
    Symbol<-rownames(object@DEData[[i]])[rownames(object@DEData[[i]]) %in% object@InputGeneSet[[i]]]
  names(FC)<-Symbol
  names(PVal)<-Symbol
  calc.object@FoldChange[[i]]<-FC
  calc.object@Pvalues[[i]]<-PVal
  UpDnCV<-vector()
  UpDnCV[calc.object@FoldChange[[i]]>0]<-"Up"
  UpDnCV[calc.object@FoldChange[[i]]<0]<-"Dn"
  calc.object@UpDnClassvec[[i]]<-as.factor(UpDnCV)
  Updf<-data.frame(Symbol=as.character(names(calc.object@FoldChange[[i]])[calc.object@FoldChange[[i]]>0]),
                   logFC=calc.object@FoldChange[[i]][calc.object@FoldChange[[i]]>0],
                   Pvalue=calc.object@Pvalues[[i]][calc.object@FoldChange[[i]]>0], stringsAsFactors=FALSE)
  Dndf<-data.frame(Symbol=names(calc.object@FoldChange[[i]])[calc.object@FoldChange[[i]]<0],
                   logFC=calc.object@FoldChange[[i]][calc.object@FoldChange[[i]]<0],
                   Pvalue=calc.object@Pvalues[[i]][calc.object@FoldChange[[i]]<0], stringsAsFactors=FALSE)
  colnames(Updf)<-c("Symbol", "logFC", "adj.P.Val")
  colnames(Dndf)<-c("Symbol", "logFC", "adj.P.Val")
  calc.object@Up[[i]]<-Updf
  calc.object@Down[[i]]<-Dndf
  names(calc.object@Up)[i]<-names(object@InputGeneSet[i])
  names(calc.object@Down)[i]<-names(object@InputGeneSet[i])
  #names(calc.object@Up[[i]]$Pvalue)<-calc.object@Up[[i]]$Symbol
  #names(calc.object@Down[[i]]$Pvalue)<-calc.object@Down[[i]]$Symbol
  if(GOopt==TRUE){
#   allGenes.Up<-object@DEData[[i]]$adj.P.Val[object@DEData[[i]]$logFC>0]
#   allGenes.Dn<-object@DEData[[i]]$adj.P.Val[object@DEData[[i]]$logFC<0]
#   names(allGenes.Up)<-rownames(object@DEData[[i]][object@DEData[[i]]$logFC>0,])
#   names(allGenes.Dn)<-rownames(object@DEData[[i]][object@DEData[[i]]$logFC<0,])
  #calc.object@GO.UP[[i]]<-runGO(allGenes=allGenes, geneSel=calc.object@Pvalues[[i]][calc.object@UpDnClassvec[[i]]=="Up"])
  #calc.object@GO.DN[[i]]<-runGO(allGenes=allGenes, geneSel=calc.object@Pvalues[[i]][calc.object@UpDnClassvec[[i]]=="Dn"])
  #calc.object@GO.UP[[i]]<-runGO(allGenes=allGenes.Up)
  #calc.object@GO.DN[[i]]<-runGO(allGenes=allGenes.Dn)
    allGenes.Up<-rep(0,length=nrow(calc.object@DEData[[i]]))
    allGenes.Up[rownames(calc.object@DEData[[i]]) %in% names(calc.object@Pvalues[[i]][calc.object@UpDnClassvec[[i]]=="Up"])]<-1
    allGenes.Up[!(rownames(calc.object@DEData[[i]]) %in% names(calc.object@Pvalues[[i]][calc.object@UpDnClassvec[[i]]=="Up"]))]<-0
    names(allGenes.Up)<-rownames(calc.object@DEData[[i]])
    allGenes.Dn<-rep(0,length=nrow(calc.object@DEData[[i]]))
    allGenes.Dn[rownames(calc.object@DEData[[i]]) %in% names(calc.object@Pvalues[[i]][calc.object@UpDnClassvec[[i]]=="Dn"])]<-1
    allGenes.Dn[!(rownames(calc.object@DEData[[i]]) %in% names(calc.object@Pvalues[[i]][calc.object@UpDnClassvec[[i]]=="Dn"]))]<-0
    names(allGenes.Dn)<-rownames(calc.object@DEData[[i]])
    calc.object@GO.BP.UP[[i]]<-runGO(allGenes=allGenes.Up, ont="BP")
    calc.object@GO.BP.DN[[i]]<-runGO(allGenes=allGenes.Dn, ont="BP")
    calc.object@GO.MF.UP[[i]]<-runGO(allGenes=allGenes.Up, ont="MF")
    calc.object@GO.MF.DN[[i]]<-runGO(allGenes=allGenes.Dn, ont="MF")
  names(calc.object@GO.BP.UP)[i]<-names(object@InputGeneSet[i])
  names(calc.object@GO.BP.DN)[i]<-names(object@InputGeneSet[i])
  names(calc.object@GO.MF.UP)[i]<-names(object@InputGeneSet[i])
  names(calc.object@GO.MF.DN)[i]<-names(object@InputGeneSet[i])
  }
  }
if(PWopt==TRUE){
  calc.object<-runDAVID(calc.object)
}
  return(calc.object)
}

topDiffGenes<-function(allScore){
  return(allScore < 0.05)
}

includeGenes<-function(allGenes){
  return(allGenes==1)
}

runGO<-function(allGenes, ont){
  if(length(allGenes)>5){
  includeGenes(allGenes)
  GOdata <- new("topGOdata",
                description = "NoRxvHC", ontology = ont,
                allGenes = allGenes, geneSel = includeGenes,
                annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

  #resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  #resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  Results<-GenTable(GOdata, topNodes = (length(GOdata@graph@nodes)/10),
                    elimKS = resultKS.elim,
                    orderBy = "elimKS")
  return(Results)}
  else{
    return("NA")
  }
}


runDAVID<-function(object){
  david<-DAVIDWebService(email="sfurlan@uw.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  if(length(object@Up)>0){
    for(i in 1:length(object@Up)){
      df<-data.frame()
      print(paste("Obtaining pathway data for genes upregulated in ", names(object@InputGeneSet)[i],sep=""))
      if(nrow(object@Up[[i]])>0){
        e2s = toTable(org.Hs.egSYMBOL)
        result<-addList(david, e2s$gene_id[e2s$symbol %in% object@Up[[i]]$Symbol],
                      idType="ENTREZ_GENE_ID",
                      listName=paste(names(object@InputGeneSet)[i], ".UP", sep=""), listType="Gene")
        setAnnotationCategories(david, c("PANTHER_PATHWAY","BIOCARTA", "REACTOME_PATHWAY","KEGG_PATHWAY"))
        Table<-getFunctionalAnnotationChart(david)
        df<-data.frame(Table@.Data, stringsAsFactors=FALSE)
        colnames(df)<-(Table@names)
        }
      if(length(df$Genes)>0){
        for(j in 1:length(df$Genes)){
          df$Genes[j]<-paste(e2s$symbol[e2s$gene_id %in% unlist(strsplit(df$Genes[j], ", "))], collapse=", ")
      }
      }
      if(nrow(object@Up[[i]])==0){
        df[1,1]<-"No genes to perform pathway analysis"
      }
      if(nrow(df)==0){
        df[1,1]<-"No pathways found"
      }
      object@PW.UP[[i]]<-df
      names(object@PW.UP)[i]<-names(object@InputGeneSet)[i]
      }
  }
  if(length(object@Down)>0){
    for(i in 1:length(object@Down)){
      df<-data.frame()
      print(paste("Obtaining pathway data for genes downregulated in ", names(object@InputGeneSet)[i],sep=""))
      if(nrow(object@Down[[i]])>0){
        e2s = toTable(org.Hs.egSYMBOL)
        result<-addList(david, e2s$gene_id[e2s$symbol %in% object@Down[[i]]$Symbol],
                    idType="ENTREZ_GENE_ID",
                    listName=paste(names(object@InputGeneSet)[i], ".DN", sep=""), listType="Gene")
        setAnnotationCategories(david, c("PANTHER_PATHWAY","BIOCARTA", "REACTOME_PATHWAY","KEGG_PATHWAY"))
        Table<-getFunctionalAnnotationChart(david)
        df<-data.frame(Table@.Data, stringsAsFactors=FALSE)
        colnames(df)<-(Table@names)
      }
      if(length(df$Genes)>0){
        for(j in 1:length(df$Genes)){
          df$Genes[j]<-paste(e2s$symbol[e2s$gene_id %in% unlist(strsplit(df$Genes[j], ", "))], collapse=", ")
        }
      }
      if(nrow(object@Down[[i]])==0){
        df[1,1]<-"No genes to perform pathway analysis"
      }
      if(nrow(df)==0){
        df[1,1]<-"No pathways found"
      }
      object@PW.DN[[i]]<-df
      names(object@PW.DN)[i]<-names(object@InputGeneSet)[i]
    }
  }
return(object)
}



writeDEAnalyzer<-function(object, directory="", cellStyle="", IPAopt=TRUE){
  if(directory!="") {
    wd<-getwd()
    newpath<-paste(wd,"/",directory, sep="")
    dir.create(newpath, showWarnings=FALSE)
  }
  files<-c("-GenesUP.xls","-GenesDN.xls", "-GenesALL.xls")
  if(object@ObjectInfo$ContainsGOData==TRUE){
    files<-c(files, "-GOTermsUP-BP.xls","-GOTermsDN-BP.xls", "-GOTermsUP-MF.xls","-GOTermsDN-MF.xls")
  }
  if(object@ObjectInfo$ContainsPWData==TRUE){
    files<-c(files, "-PathwaysUP.xls","-PathwaysDN.xls")
  }
  #data.ind<-c("-GenesUP.xls","-GenesDN.xls", "-GOTermsUP.xls","-GOTermsDN.xls", "-PathwaysUP.xls","-PathwaysDN.xls") %in% files
  prefix<-format(Sys.Date(), format="%Y%m%d")
  if(directory!="") {
    location<-paste(directory, "/",sep="")
    filename<-prefix
  }
  else {
    filename<-prefix
    location<-""
  }
  if(IPAopt==TRUE){
    if(directory!="") {
      newpath<-paste(newpath,"/IPA", sep="")
      dir.create(newpath, showWarnings=FALSE)
      IPAlocation<-paste(location, "IPA/", sep="")
    }
    else{IPAlocation<-"IPA/"}
  }
  wb<-list()
  for(j in 1:length(files)){
    wb[[j]]<-loadWorkbook(filename=paste(location, filename, files[j], sep=""), create=TRUE)
  }
  for(i in 1:length(object@InputGeneSet)){
    dflist<-list(Genesup.df<-object@Up[[i]],
         Genesdn.df<-object@Down[[i]], Genesall.df<-rbind(object@Up[[i]], object@Down[[i]]))
        if(object@ObjectInfo$ContainsGOData==TRUE){
          dflist<-c(dflist,
                    list(dflistGOup.df<-object@GO.BP.UP[[i]]), list(GOdn.df<-object@GO.BP.DN[[i]]),
                    list(dflistGOup.df<-object@GO.MF.UP[[i]]), list(GOdn.df<-object@GO.MF.DN[[i]]))
        }
        if(object@ObjectInfo$ContainsPWData==TRUE){
          dflist<-c(dflist, list(dflistPWup.df<-object@PW.UP[[i]]), list(PWdn.df<-object@PW.DN[[i]]))
        }
    if(IPAopt==TRUE){
      IPAout<-cbind(rownames(object@DEData[[i]]), object@DEData[[i]][,c(1,5)])
      colnames(IPAout)<-c("Symbol", "Log.Ratio", "Adj.P.Val")
      write.table(IPAout, file=paste(IPAlocation, filename, "-", names(object@InputGeneSet)[i], ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
    }
    for(j in 1:length(files)){
      createSheet(wb[[j]], name = names(object@InputGeneSet)[i])
      writeWorksheet(wb[[j]], dflist[[j]], sheet=names(object@InputGeneSet)[i])
#       if(class(cellstyle)!="character"){
#       createCellStyle(wb[[j]], "GeneData")
#       setCellStyle(wb[[j]], sheet =names(object@InputGeneSet)[i], row =1 , col =1 , cellstyle =cellStyle)
#       }
      saveWorkbook(wb[[j]])
    }
  }
}

quickXL<-function(filename=paste(format(Sys.Date(), format="%Y_%m_%d_"), format(Sys.time(), format="%H_%M"), ".xls", sep=""), obj, sheetname="1"){
  if(class(obj)!="data.frame"){stop("Only data frames supported")}
     wb<-loadWorkbook(filename=filename, create=TRUE)
    createSheet(wb, name = sheetname)
  obj2<-cbind(rownames(obj),obj)
  colnames(obj2)<-c("Symbol", colnames(obj))
      writeWorksheet(wb, obj2, sheet=sheetname)
     saveWorkbook(wb)
}


quickXL.list<-function(filename=NULL, obj, header=FALSE){
  if(class(obj)!="list"){stop("Only lists supported")}
  if(is.null(filename)){stop("Invalid filename")}
  wb<-loadWorkbook(filename=filename, create=TRUE)
  for(i in names(obj)){
  createSheet(wb, name=i)
  dat<-obj[[i]]
  writeWorksheet(wb, dat, sheet=i, header=header)
  }
  saveWorkbook(wb)
}



readautoLIMMA<-function(limmalist){
toprint<-paste(1:length(as.character(limmalist[[2]][1,])), as.character(limmalist[[2]][1,]), sep="-")
return(toprint)
}

autoDAVID<-function(vector, input="symbol"){
  if(symbol %in% c("symbol", "ENS")){break("Input Not Correct")}
  library(org.Hs.eg.db)
  david<-RDAVIDWebService::DAVIDWebService(email="sfurlan@uw.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  if(input=="symbol"){
  for(i in 1:length(vector)){
    e2s = toTable(org.Hs.egSYMBOL)
    result<-addList(david, e2s$gene_id[e2s$symbol %in% vector],
                    idType="ENTREZ_GENE_ID",
                    listName="autoDavid", listType="Gene")
    found<-length(which(e2s$gene_id %in% vector))
    print(paste("Found ", found, " Annotated Genes of ", length(vector), " Submitted", sep=""))
    }
  }
  if(input=="ENS"){
    e2s <- toTable(org.Hs.egENSEMBLTRANS)
    result<-addList(david, e2s$gene_id[e2s$trans_id %in% vector],
                    idType="ENTREZ_GENE_ID",
                    listName="autoDavid", listType="Gene")
    found<-length(which(e2s$trans_id %in% vector))
    print(paste("Found ", found, " Annotated Genes of ", length(vector), " Submitted", sep=""))
  }
    setAnnotationCategories(david, c("PANTHER_PATHWAY","BIOCARTA", "REACTOME_PATHWAY","KEGG_PATHWAY"))
    Table<-getFunctionalAnnotationChart(david)
    df<-data.frame(Table@.Data, stringsAsFactors=FALSE)
    colnames(df)<-(Table@names)
    if(symbol=="symbol"){
    if(length(df$Genes)>0)
      {
        for(j in 1:length(df$Genes)){
          df$Genes[j]<-paste(e2s$symbol[e2s$gene_id %in% unlist(strsplit(df$Genes[j], ", "))], collapse=", ")
        }
      }
    }
  return(df)
}

GroupSelection<-function(classvec, SelectionIndex, data, alldata, cols, method="separate"){
  selected<<-levels(classvec)[SelectionIndex]
  data.sel<<-data[,classvec %in% selected]
  classvec.sel<<-as.factor(as.character(classvec[classvec %in% selected]))
  alldata.sel<<-alldata[,classvec %in% selected]
  colnames(alldata.sel)<-make.unique(as.character(classvec.sel))
  colmatch<-cols$cols[match(selected,cols$group)]
  names(colmatch)<-cols$group[match(selected,cols$group)]
  #pie(rep(1,length(colmatch)), col=colmatch, labels=names(colmatch))
  gcolor<<-colmatch[order(names(colmatch))]
  pie(rep(1,length(gcolor)), col=gcolor, labels=names(gcolor))
  tcolor<<-gcolor[names(gcolor)!="Healthy.Control"]
  gcolors<<-gcolor[classvec.sel]
  tcolors<<-tcolor[classvec[classvec.sel!="Healthy.Control"]]
  limma.master<<-autoLIMMA(data.sel, classvec.sel, element1=seq(1, length(levels(classvec.sel))), element2=seq(1, length(levels(classvec.sel))), pvalue.thresh=0.05, lfc.thresh=0.5, adjust.method="BH", method=method)
  limma.alldata<<-autoLIMMA(alldata.sel, classvec.sel, element1=seq(1, length(levels(classvec.sel))), element2=seq(1, length(levels(classvec.sel))), pvalue.thresh=0.05, lfc.thresh=0.5, adjust.method="BH", method=method)
#   return(invisible(limma.master))
#   return(invisible(limma.alldata))
#   return(invisible(gcolors))
#   return(invisible(tcolors))
#   return(invisible(tcolor))
#   return(invisible(gcolor))
#   return(invisible(data.sel))
#   return(invisible(alldata.sel))
#   return(invisible(classvec.sel))
}

GroupSelectionNoLimma<-function(classvec, SelectionIndex, data, alldata, cols){
  selected<<-levels(classvec)[SelectionIndex]
  data.sel<<-data[,classvec %in% selected]
  classvec.sel<<-as.factor(as.character(classvec[classvec %in% selected]))
  alldata.sel<<-alldata[,classvec %in% selected]
  colnames(alldata.sel)<-make.unique(as.character(classvec.sel))
  colmatch<-cols$cols[match(selected,cols$group)]
  names(colmatch)<-cols$group[match(selected,cols$group)]
  #pie(rep(1,length(colmatch)), col=colmatch, labels=names(colmatch))
  gcolor<<-colmatch[order(names(colmatch))]
  pie(rep(1,length(gcolor)), col=gcolor, labels=names(gcolor))
  tcolor<<-gcolor[names(gcolor)!="Healthy.Control"]
  gcolors<<-gcolor[classvec.sel]
  tcolors<<-tcolor[classvec[classvec.sel!="Healthy.Control"]]
#   limma.master<<-autoLIMMA(data.sel, classvec.sel, element1=seq(1, length(levels(classvec.sel))), element2=seq(1, length(levels(classvec.sel))), pvalue.thresh=0.05, lfc.thresh=0.5, adjust.method="BH")
#   limma.alldata<<-autoLIMMA(alldata.sel, classvec.sel, element1=seq(1, length(levels(classvec.sel))), element2=seq(1, length(levels(classvec.sel))), pvalue.thresh=0.05, lfc.thresh=0.5, adjust.method="BH")
#   #   return(invisible(limma.master))
  #   return(invisible(limma.alldata))
  #   return(invisible(gcolors))
  #   return(invisible(tcolors))
  #   return(invisible(tcolor))
  #   return(invisible(gcolor))
  #   return(invisible(data.sel))
  #   return(invisible(alldata.sel))
  #   return(invisible(classvec.sel))
}

GroupSelectionLocal<-function(suff=".new", classvec, SelectionIndex, data, alldata, cols, method="separate"){
  selected<-levels(classvec)[SelectionIndex]
  data.sel<-data[,classvec %in% selected]
  classvec.sel<-as.factor(as.character(classvec[classvec %in% selected]))
  alldata.sel<-alldata[,classvec %in% selected]
  assign(paste("selected", suff, sep=""), selected, envir=.GlobalEnv)
   assign(paste("data.sel", suff, sep=""), data.sel, envir=.GlobalEnv)
   assign(paste("classvec.sel", suff, sep=""), classvec.sel, envir=.GlobalEnv)
   assign(paste("alldata.sel", suff, sep=""), alldata.sel, envir=.GlobalEnv)
   colnames(alldata.sel)<-make.unique(as.character(classvec.sel))
  colmatch<-cols$cols[match(selected,cols$group)]
  names(colmatch)<-cols$group[match(selected,cols$group)]
  gcolor<-colmatch[order(names(colmatch))]
  pie(rep(1,length(gcolor)), col=gcolor, labels=names(gcolor))
  tcolor<-gcolor[names(gcolor)!="Healthy.Control"]
  gcolors<-gcolor[classvec.sel]
  tcolors<-tcolor[classvec[classvec.sel!="Healthy.Control"]]
  limma.master<-autoLIMMA(data.sel, classvec.sel, element1=seq(1, length(levels(classvec.sel))), element2=seq(1, length(levels(classvec.sel))), pvalue.thresh=0.05, lfc.thresh=0.5, adjust.method="BH", method=method)
  limma.alldata<-autoLIMMA(alldata.sel, classvec.sel, element1=seq(1, length(levels(classvec.sel))), element2=seq(1, length(levels(classvec.sel))), pvalue.thresh=0.05, lfc.thresh=0.5, adjust.method="BH", method=method)
  assign(paste("gcolor", suff, sep=""), gcolor, envir=.GlobalEnv)
  assign(paste("tcolor", suff, sep=""), tcolor, envir=.GlobalEnv)
  assign(paste("gcolors", suff, sep=""), gcolors, envir=.GlobalEnv)
  assign(paste("tcolors", suff, sep=""), tcolors, envir=.GlobalEnv)
  assign(paste("limma.master", suff, sep=""), limma.master, envir=.GlobalEnv)
  assign(paste("limma.alldata", suff, sep=""), limma.alldata, envir=.GlobalEnv)
}
#debug(AssignColors)
AssignColors<-function(classvec, cols){
  selected<-levels(classvec)
  colmatch<-lapply(cols, "[", match(selected,cols$group))
  #names(colmatch)<-cols[match(selected,cols$group)]
  #pie(rep(1,length(colmatch)), col=colmatch, labels=names(colmatch))
  #gcolor<<-colmatch[order(names(colmatch))]
  pie(rep(1,length(selected)), col=colmatch$cols, labels=colmatch$group)
  #tcolor<<-gcolor[names(gcolor)!="Healthy.Control"]
  #gcolors<<-gcolor[classvec]
  #tcolors<<-tcolor[classvec[classvec!="Healthy.Control"]]
  return(colmatch)
}
