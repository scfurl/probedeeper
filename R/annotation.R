# vec<-  c("A_23_P111206", "A_23_P126540", "A_23_P145965", "A_23_P146233", "A_23_P17065",  "A_23_P206120", "A_23_P259071",
# "A_23_P310460", "A_23_P321307", "A_23_P326963", "A_23_P362719", "A_23_P56494",  "A_23_P79398",  "A_23_P85082",
# "A_23_P8640",   "A_23_P86470" , "A_23_P92672",  "A_23_P99442" , "A_24_P230116", "A_24_P251866", "A_24_P257224",
# "A_24_P3016" ,  "A_24_P63019" , "A_24_P686247", "A_24_P917668", "A_24_P96505" , "A_32_P113190", "A_32_P142521",
# "A_32_P170879", "A_32_P17743")
# file<-"~/Dropbox (Kean Lab)/AWS/Scott/Rproj/SampleData/GSE23924/GSE23924_RAW/GPL6480.soft"
# file.exists(file)


annotateProbes<-function(vec, file=NULL, filetype="csv"){
  if(length(file) == 0) { break("No File input") }
  if(class(vec) != "character") { break("Not correct input") }
  if(filetype == "csv") { input<-read.csv(file) }
  if(filetype == "soft") { input<-read.csv(file) }
  input.soft<-GEOquery::getGEO(filename=file)
  annt<-loadGeoFile(file)
  firstpass<-annt[vec,"GENE_SYMBOL"]
  output<-vec
  output[firstpass!=""]<-firstpass[firstpass!=""]
  #output[firstpass==""]<-annt[output[firstpass==""],"REFSEQ"]
  return(output)
}




masterannotate<-function(obj1, gene.multiples=FALSE, delete.NA=TRUE, method="maxRowVariance"){
  require(WGCNA)
  ###gene.multiples = t/f allow multiple gene annotations - will make unique symbol names
  ###delete.NA = delete NA annotations; if FALSE will annotate with probe names
  ### This function removes unannotated genes ####
  ### This function also calculates median expression for genes that have multiple probes ####
  annotation.file<-read.csv(paste(rootpath,"/Bioinformatics Resources/Rhesus Annotation/MasterAnnotation.csv", sep=""), header= TRUE, colClasses='character')
  symbol<-data.frame(Symbol=annotation.file$Symbol, row.names=annotation.file$ID, stringsAsFactors=FALSE)
  #colnames(symbol)<-"Symbol"
  ags<-as.data.frame(symbol[row.names(obj1),], stringsAsFactors=FALSE)
  colnames(ags)<-"Symbol"
  obj1.an<-as.data.frame(obj1, stringsAsFactors=FALSE)
  obj1.an$Symbol<-ags$Symbol
  obj2.an<-obj1.an[complete.cases(obj1.an),]
  obj2.pr<-obj1.an[is.na(obj1.an$Symbol==TRUE),]
  ###For no selection and all data preserved###
  if(delete.NA==FALSE && gene.multiples ==TRUE){
    rownames(obj2.an)<-make.unique(obj2.an$Symbol, sep = "#")
    obj2.an$Symbol<-NULL
    obj2.pr$Symbol<-NULL
    final<-rbind(obj2.an, obj2.pr)
    return(final)}
  ###For no selection and eliminate NA###
  if(delete.NA==TRUE && gene.multiples ==TRUE){
    rownames(obj2.an)<-make.unique(obj2.an$Symbol, sep = "#")
    obj2.an$Symbol<-NULL
    return(obj2.an)}
  ####For selection and eliminated NA####
  if(delete.NA==TRUE && gene.multiples ==FALSE){
    sel.obj.an<-obj2.an
    Symbol<-sel.obj.an$Symbol
    sel.obj.an$Symbol<-NULL
  }
  if(delete.NA==FALSE && gene.multiples ==FALSE){
    obj2.pr$Symbol<-rownames(obj2.pr)
    sel.obj.an<-rbind(obj2.an, obj2.pr)
    Symbol<-sel.obj.an$Symbol
    sel.obj.an$Symbol<-NULL
  }
  ###SELECTION FUNCTIONS###
  collapse.object=collapseRows(datET=sel.obj.an, rowGroup=Symbol, rowID=rownames(sel.obj.an), method=method)
  return(collapse.object$datETcollapsed)
}

annotatePerFile<-function(obj1, file=NULL, gene.multiples=FALSE, delete.NA=TRUE, method="maxRowVariance"){
  require(WGCNA)
  ###gene.multiples = t/f allow multiple gene annotations - will make unique symbol names
  ###delete.NA = delete NA annotations; if FALSE will annotate with probe names
  ### This function removes unannotated genes ####
  ### This function also calculates median expression for genes that have multiple probes ####
  annotation.file<-read.csv(file, header= TRUE, colClasses='character')
  symbol<-data.frame(Symbol=annotation.file$Symbol, row.names=annotation.file$ID, stringsAsFactors=FALSE)
  #colnames(symbol)<-"Symbol"
  ags<-as.data.frame(symbol[row.names(obj1),], stringsAsFactors=FALSE)
  colnames(ags)<-"Symbol"
  obj1.an<-as.data.frame(obj1, stringsAsFactors=FALSE)
  obj1.an$Symbol<-ags$Symbol
  obj2.an<-obj1.an[complete.cases(obj1.an),]
  obj2.pr<-obj1.an[is.na(obj1.an$Symbol==TRUE),]
  ###For no selection and all data preserved###
  if(delete.NA==FALSE && gene.multiples ==TRUE){
    rownames(obj2.an)<-make.unique(obj2.an$Symbol, sep = "#")
    obj2.an$Symbol<-NULL
    obj2.pr$Symbol<-NULL
    final<-rbind(obj2.an, obj2.pr)
    return(final)}
  ###For no selection and eliminate NA###
  if(delete.NA==TRUE && gene.multiples ==TRUE){
    rownames(obj2.an)<-make.unique(obj2.an$Symbol, sep = "#")
    obj2.an$Symbol<-NULL
    return(obj2.an)}
  ####For selection and eliminated NA####
  if(delete.NA==TRUE && gene.multiples ==FALSE){
    sel.obj.an<-obj2.an
    Symbol<-sel.obj.an$Symbol
    sel.obj.an$Symbol<-NULL
  }
  if(delete.NA==FALSE && gene.multiples ==FALSE){
    obj2.pr$Symbol<-rownames(obj2.pr)
    sel.obj.an<-rbind(obj2.an, obj2.pr)
    Symbol<-sel.obj.an$Symbol
    sel.obj.an$Symbol<-NULL
  }
  ###SELECTION FUNCTIONS###
  collapse.object=collapseRows(datET=sel.obj.an, rowGroup=Symbol, rowID=rownames(sel.obj.an), method=method)
  return(collapse.object$datETcollapsed)
}


loadGeoFile <- function(geoFilename) {
  temp <- readLines(geoFilename)    # Load the file
  temp <- temp[grep("\t", temp)]    # Keep only lines with tabs
  temp <- gsub("\t$", "\tNA", temp) # Deal with NA
  temp <- strsplit(temp, "\t")      # Split the strings at each tab
  temp <- t(sapply(temp, unlist))   # Turn each line into a vector, transpose
  colnames(temp) <- temp[1, ]
  rownames(temp) <- temp[ ,1]
  #Remove the row/col names from the data, and return it.
  #Note that all the entries are strings/characters, not numeric!
  temp[-1,-1]
}
