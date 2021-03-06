
LimmaObjCreate<-function(eset,
                         ColObj,
                         element1=NULL,
                         element2=NULL, pvalue.thresh=0.05, lfc.thresh=1, adjust.method="fdr", method="separate", printdata=FALSE, norm.method="unspecified"){
  ####Auto-LIMMA####
  #inputs data.sel, classvec.sel, element1,(2), pvalue.thresh, lfc.thresh
  classvec.sel<-ColObj@classvec
  if(is.null(element1) && is.null(element2)){
    element1<-1:length(levels(classvec.sel))
    element2<-element1}
  data.sel<-exprs(eset)
  limmaDE.output<-list()
  mydesign <- model.matrix(~0 + classvec.sel)
  LETTERS2600 <- c(sapply(LETTERS, function(x) paste0(x, 1:100)))
  colnames(mydesign) <- LETTERS2600[seq(from=1, to=length(levels(classvec.sel)))]
  names(colnames(mydesign))<-levels(classvec.sel)
  fit <- limma::lmFit(data.sel, mydesign)
  contrast.string<-contrast.eq(colnames(mydesign), element1=element1, element2=element2)
#   args<-as.list(c(contrast.string, "mydesign"))
#   names(args)[length(args)]<-"levels"
  contrast.matrix<-limma::makeContrasts(contrasts=contrast.string, levels=mydesign)
  space1<-unlist(lapply(gregexpr(" ", contrast.string), "[[", 1))
  space2<-unlist(lapply(gregexpr(" ", contrast.string), "[[", 2))
  space3<-unlist(lapply(gregexpr(" ", contrast.string), "[[", 3))
  space4<-unlist(lapply(gregexpr(" ", contrast.string), "[[", 4))
  colnames(contrast.matrix)<-substring(contrast.string,1,space1-1)
  genofits <- limma::contrasts.fit(fit, contrast.matrix)
  geno_ebFit <- limma::eBayes(genofits)
  results<-limma::decideTests(geno_ebFit, method)
  contrast.list<-as.list(substring(contrast.string,1,space1-1))
  #contrast.list<-as.list(substring(contrast.string,1,2))
  names(contrast.list)<-as.list(substring(contrast.string,1,space1-1))
  contrast.mat<-as.data.frame(mapply(append, contrast.list, substring(contrast.string,space2+1, space3-1)), stringsAsFactors=FALSE)
  contrast.mat<-rbind(contrast.mat, substring(contrast.string,space4+1))
  colnames(contrast.mat)<-paste("Contrast", as.character(1:ncol(contrast.mat)), sep="")
  tmp<-paste(names(colnames(mydesign))[match(as.character(contrast.mat[2,]), colnames(mydesign))],
             names(colnames(mydesign))[match(contrast.mat[3,], colnames(mydesign))], sep="v")
  contrast.mat<-rbind(contrast.mat, tmp)
  contrast.list<-list(coef=contrast.mat[1,], meaning=contrast.mat[4,])
  contrast.df<-data.frame(
    label=colnames(contrast.mat),
    coef=as.character(contrast.mat[1,]),
    meaning=as.character(contrast.mat[4,]),
    comp1=names(colnames(mydesign))[match(contrast.mat[2,], colnames(mydesign))],
    comp2=names(colnames(mydesign))[match(contrast.mat[3,], colnames(mydesign))],
    stringsAsFactors=FALSE)
  limmaDE.output<-contrast.list
  LimmaObj<-new('LimmaObj')
  tmp.list<-list()
  tmp.list2<-list()
  for(i in 1:length(contrast.list[[1]])){
    tmp.list2[[i]]<-limma::topTable(geno_ebFit,coef=as.character(limmaDE.output[[1]][i]),n=60000,adjust.method=adjust.method, sort.by="logFC", p.value=pvalue.thresh,lfc=lfc.thresh)
    names(tmp.list2)[i]<-paste(as.character(limmaDE.output[[1]][i]), "-", as.character(limmaDE.output[[2]][i]), sep="")
    #if(os=="Windows"){rownames(tmp.list2[[i]])<-tmp.list2[[i]]$ID}
    tmp.list[[i]]<-limma::topTable(geno_ebFit,coef=as.character(limmaDE.output[[1]][i]),n=60000,adjust.method=adjust.method, sort.by="logFC")
    names(tmp.list)[i]<-paste(as.character(limmaDE.output[[1]][i]), "-", as.character(limmaDE.output[[2]][i]), sep="")
    #if(os=="Windows"){rownames(tmp.list[[i]])<-tmp.list[[i]]$ID}
  }
  LimmaObj@DEGenes<-tmp.list2
  LimmaObj@AllGenes<-tmp.list
  #   names(limmaDE.output)[4]<-"AllLimmaDEData"
  #   names(limmaDE.output)[3]<-"SigLimmaDEData"

  LimmaObj@Contrasts<-contrast.df
  #   limmaDE.output[5]<-list(Comparisons.df=data.frame(Comp1=comp1, Comp2=comp2))
  LimmaObj@Inputs<-list(p.adjust=adjust.method, p.thresh=pvalue.thresh, lfc.thresh=lfc.thresh, norm.method=norm.method)
  #   if(printdata==TRUE){
  #     print(names(limmaDE.output[[3]]))
  #   }
  return(LimmaObj)
}

ExtractLimmaObj<-function(gene, LimmaObj){
  rownamesdf<-LimmaObj@Contrasts$meaning
  sigLevel<-vector()
  logFC<-vector()
  PVal<-vector()
  for(i in 1:length(rownamesdf)){
    logFC[i]<-LimmaObj@AllGenes[[i]][gene,]$logFC
    PVal[i]<-LimmaObj@AllGenes[[i]][gene,]$adj.P.Val
    sigLevel[i]<-SigLevel(PVal[i])
  }
  df<-data.frame(logFC=logFC, adj.P.Val=PVal, sigLevel=sigLevel, row.names=rownamesdf)
  return(df)
}


jgc <- function()
{
  .jcall("java/lang/System", method = "gc")
}


LimmaObj2XL<-function(object,
                      directory=get.wd(),
                      prefix="LimmaObj",
                      type="DEGenes",
                      index=NULL,
                      method="openxlsx"){
  if(class(object)!="LimmaObj"){stop("Input does not appear to be a LimmaObj")}
  files<-paste(prefix, "-", c("GenesUP.xls","GenesDN.xls", "GenesALL.xls"), sep="")
  wb<-list()
  if(is.null(index)){
    index<-1:length(object@Contrasts$meaning)
  }
  df<-data.frame(str=object@Contrasts$meaning[index], index=index)
  if(method=="XLConnect"){
    for(j in 1:length(files)){
      wb[[j]]<-XLConnect::loadWorkbook(filename=file.path(directory, files[j]), create=TRUE)
    }
  }else{
    for(j in 1:length(files)){
      wb[[j]]<-createWorkbook()
    }
  }
  if(type=="AllGenes"){
    for(i in df$index){
      dflist<-list(Genesup.df<-object@AllGenes[[i]][object@AllGenes[[i]]$logFC>0,],
                   Genesdn.df<-object@AllGenes[[i]][object@AllGenes[[i]]$logFC<0,], Genesall.df<-object@AllGenes[[i]])
      if(method=="XLConnect"){
        for(j in 1:length(files)){
          jgc()
          message("Creating sheet ", object@Contrasts$meaning[i])
          createSheet(wb[[j]], name = object@Contrasts$meaning[i])
          writeWorksheet(wb[[j]], dflist[[j]], sheet=object@Contrasts$meaning[i], rownames="Symbol")
          #       if(class(cellstyle)!="character"){
          #       createCellStyle(wb[[j]], "GeneData")
          #       setCellStyle(wb[[j]], sheet =names(object@InputGeneSet)[i], row =1 , col =1 , cellstyle =cellStyle)
          #       }
          saveWorkbook(wb[[j]])
          xlcFreeMemory()
          xlcMemoryReport()
        }
      }else{
        for(j in 1:length(files)){
          message("Creating sheet ", object@Contrasts$meaning[i])
          addWorksheet(wb[[j]], sheetName = object@Contrasts$meaning[i])
          writeData(wb = wb[[j]], sheet = object@Contrasts$meaning[i], x = dflist[[j]], borders = "n", rowNames = TRUE)
          saveWorkbook(wb[[j]], file.path(directory, files[j]), overwrite = TRUE)
        }
      }
    }
  }
  if(type=="DEGenes"){
    for(i in df$index){
      dflist<-list(Genesup.df<-object@DEGenes[[i]][object@DEGenes[[i]]$logFC>0,],
                   Genesdn.df<-object@DEGenes[[i]][object@DEGenes[[i]]$logFC<0,], Genesall.df<-object@DEGenes[[i]])
      if(method=="XLConnect"){
        for(j in 1:length(files)){
          createSheet(wb[[j]], name = object@Contrasts$meaning[i])
          writeWorksheet(wb[[j]], dflist[[j]], sheet=object@Contrasts$meaning[i], rownames="Symbol")
          #       if(class(cellstyle)!="character"){
          #       createCellStyle(wb[[j]], "GeneData")
          #       setCellStyle(wb[[j]], sheet =names(object@InputGeneSet)[i], row =1 , col =1 , cellstyle =cellStyle)
          #       }
          saveWorkbook(wb[[j]])
        }
      }else{
        for(j in 1:length(files)){
          message("Creating sheet ", object@Contrasts$meaning[i])
          addWorksheet(wb[[j]], sheetName = object@Contrasts$meaning[i])
          writeData(wb = wb[[j]], sheet = object@Contrasts$meaning[i], x = dflist[[j]], borders = "n", rowNames = TRUE)
          saveWorkbook(wb[[j]], file.path(directory, files[j]), overwrite = TRUE)
        }
      }
    }
  }
}


# makeContrasts.SF<-function(contrast.string) {
#   contrast.string=c(contrast.string, "levels = mydesign")
#   contrasts<-do.call(limma::makeContrasts, as.list(contrast.string), quote=TRUE)
#   return(contrasts)
# }
