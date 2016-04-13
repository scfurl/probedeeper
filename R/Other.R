#####OTHER FUNCTIONS########

synantag<-function(pred, act){
  if(length(pred)!=length(act)){stop("data is not of same length")}
  relationship<-character()
for(i in 1: length(pred)){
  
  if(pred[i]<0 && act[i]<0 && pred[i]>act[i]){relationship[i]<-"synergistic"}
  if(pred[i]>0 && act[i]>0 && pred[i]<act[i]){relationship[i]<-"synergistic"}
  if(isNA(relationship[i])==TRUE){relationship[i]<-"antagonistic"}}
return(relationship)}
  

###CONVERT DATA TO ENTREZ GENE###

entrezdf<-function(df){
  e2s = toTable(org.Hs.egSYMBOL)
  df.sliced<-df[rownames(df) %in% e2s$symbol,]
  rownames(df.sliced)<-e2s$gene_id[match(rownames(df.sliced), e2s$symbol)]
  return(df.sliced)
}