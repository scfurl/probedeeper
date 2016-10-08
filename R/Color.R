# file<-"~/Dropbox (Kean Lab)/AWS/Scott/Bioinformatics Resources/Color/CBSafe15.csv"
# file.exists(file)
# read.csv(file, header=TRUE)
readColorFile<-function(file, type="csv", assign=FALSE){
  if(type=="csv"){
    if(!file.exists(file)){stop("File doesn't appear to exist")}
    else{
      colors<-vector()
      data<-read.csv(file, header=TRUE)
      for(i in 1:nrow(data)){
        colors<-c(colors, rgb(data[i,2], data[i,3], data[i,4], maxColorValue=255))
      }
      names(colors)<-data[,1]
    }
    if(assign==TRUE){
      colorsubset<-colors[data$Assign!=""]
      names(colorsubset)<-data$Assign[data$Assign!=""]
      return(colorsubset)
    }
  }
  else{
    stop("Only CSV currently supported")
  }
  return(colors)
}

# colors<-readColorFile(file)
# par(mfrow=c(1,2))
# pie(rep(1, length(colors)), col=colors, labels=names(colors))
# pie(rep(1, length(colors)), col=colors, labels=1:15)

gsea.col<- c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000")
