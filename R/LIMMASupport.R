#LIMMA SUPPORT

fixrownames<-function(df, name)
{
  #This functions creates a column called "name" out of rownames and puts it in the front of a df
  output<-cbind(rownames(df), df,stringsAsFactors=FALSE)
  colnames(output)<-c(name,colnames(df))
  return(output)
}


contrast.eq<-function(vector, element1=seq(1,length(vector)), element2){
  #This function provides all possible contrants between two vectors
  grid<-expand.grid(vector[element1], vector[element2], stringsAsFactors=FALSE)
  grid[grid$Var1==grid$Var2,]<-NA
  grid<-grid[complete.cases(grid),]
  vector.out<-vector(length=nrow(grid))
  for(i in 1:nrow(grid)){
    vector.out[i]<-paste(grid[i,1], grid[i,2], " = ",grid[i,1], " - " ,grid[i,2], sep="")
  }
  return(vector.out)
}

makeContrasts.SF<-function(...) {
  limma::makeContrasts(...,levels = mydesign)
}

autoLIMMA<-function(data.sel, classvec.sel, element1, element2, pvalue.thresh=0.05, lfc.thresh=1, adjust.method="fdr", method="separate", printdata=FALSE){
####Auto-LIMMA####
#inputs data.sel, classvec.sel, element1,(2), pvalue.thresh, lfc.thresh
limmaDE.output<-list()
mydesign <- model.matrix(~0 + classvec.sel)
colnames(mydesign) <- LETTERS[seq(from=1, to=length(levels(classvec.sel)))]
names(colnames(mydesign))<-levels(classvec.sel)
fit <- limma::lmFit(data.sel, mydesign)
mydesign<<-mydesign
contrast.string<-contrast.eq(colnames(mydesign), element1=element1, element2=element2)
contrast.matrix<-do.call(makeContrasts.SF, as.list(contrast.string))
colnames(contrast.matrix)<-substring(contrast.string,1,2)
genofits <- limma::contrasts.fit(fit, contrast.matrix)
geno_ebFit <<- limma::eBayes(genofits)
results<<-limma::decideTests(geno_ebFit, method)
contrast.list<-as.list(substring(contrast.string,1,2))
names(contrast.list)<-as.list(substring(contrast.string,1,2))
contrast.mat<-as.data.frame(mapply(append, contrast.list, substring(contrast.list,1,1)), stringsAsFactors=FALSE)
contrast.mat<-rbind(contrast.mat, substring(contrast.mat[1,],2,2))
tmp<-paste(names(colnames(mydesign))[match(contrast.mat[2,], colnames(mydesign))],
           names(colnames(mydesign))[match(contrast.mat[3,], colnames(mydesign))], sep="v")
contrast.mat<-rbind(contrast.mat, tmp)
contrast.list<-list(coef=contrast.mat[1,], meaning=contrast.mat[4,])
limmaDE.output<-contrast.list
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
limmaDE.output[[3]]<-tmp.list2
limmaDE.output[[4]]<-tmp.list
names(limmaDE.output)[4]<-"AllLimmaDEData"
names(limmaDE.output)[3]<-"SigLimmaDEData"
comp1<-names(colnames(mydesign))[match(contrast.mat[2,], colnames(mydesign))]
comp2<-names(colnames(mydesign))[match(contrast.mat[3,], colnames(mydesign))]
limmaDE.output[5]<-list(Comparisons.df=data.frame(Comp1=comp1, Comp2=comp2))
limmaDE.output[[6]]<-list(p.adjust=adjust.method, p.thresh=pvalue.thresh, lfc.thresh=lfc.thresh)
if(printdata==TRUE){
  print(names(limmaDE.output[[3]]))
}
return(limmaDE.output)
}


SFVenn<-function(limma.out, index.sel, colors=RColorBrewer::brewer.pal(12, "Paired"), classvec.sel=classvec.sel, plottype = "ChowRuskey"){
  #This function takes output from autoLimma and makes Venn diagrams of the various contrasts
  #limma.out= output from autoLIMMA
  #index.sel= vector of indices to plot Venn
  #colors= vector of colors for index.sel
  VList=list()
  for(i in 1:length(index.sel)){
    VList[[i]]<- rownames(limma.out[[3]][[index.sel[i]]])
  }
  names(VList)<-limma.out[[2]][index.sel]
  Vstem <- Venn(VList)
  Tstem <- compute.Venn(Vstem)
  gp <- VennThemes(Tstem, colourAlgorithm = "sequential", increasingLineWidth=FALSE)
  gps <- gp[["Set"]]
  nSets <- length(gps)
  venncolors<-colors[which(LETTERS[seq(from=1, to=length(colors))] %in% substr(names(limma.out[[3]]),1,1)[index.sel])]
  if(length(venncolors)<nSets){
    venncolors<-c(venncolors, RColorBrewer::brewer.pal(12, "Paired")[1:(nSets-length(venncolors))])
  }
  make.unique(venncolors)
  for (ix in 1:nSets) {
    gps[[ix]]$col<-venncolors[ix]
  }
  gp[["Set"]] <- gps
  for (ix in 1:nSets) {
    gp$SetText[[ix]]$col<-venncolors[ix]
    gp$SetText[[ix]]$fontsize<-16
  }
  plot(Vstem, type=plottype, show = list(SetLabels = TRUE), gp=gp)
  for(i in 1:length(VList)){
    print(paste(names(VList[i]), length(VList[[i]]), sep="-"))
  }
  return(Vstem)
  }


ExtractLIMMA<-function(object, gene){
  rownamesdf<-as.character(apply(object[[2]],1, as.character))
  sigLevel<-vector()
  logFC<-vector()
  PVal<-vector()
  for(i in 1:length(object[[4]])){
    logFC[i]<-object[[4]][[i]][gene,]$logFC
    PVal[i]<-object[[4]][[i]][gene,]$adj.P.Val
    sigLevel[i]<-SigLevel(PVal[i])
  }
  df<-data.frame(logFC=logFC, adj.P.Val=PVal, sigLevel=sigLevel, row.names=rownamesdf)
  return(df)
}
