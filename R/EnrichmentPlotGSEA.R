EnrichmentPlotLite<-function(filelist,title, plotcols=rep("black", length(filelist)), stats){
  Rep<-read.delim(filelist[i], header=TRUE, stringsAsFactors=FALSE)
  df<-data.frame()
df<-cbind.data.frame(as.numeric(Rep$RES), as.numeric(Rep$LIST.LOC), as.character(Rep$CORE_ENRICHMENT))
colnames(df)<-c("RES", "LOC", "ENRICHMENT")
df$max<-df$LOC[which(df$RES == max(df$RES))][1]
df$min<-df$LOC[which(df$RES == min(df$RES))][1]

lRES[i]<-df$RES[which(df$RES == min(df$RES))][1]
uRES[i]<-df$RES[which(df$RES == max(df$RES))][1]
if(abs(lRES[i])>(abs(uRES[i]))){df$dotted<-df$min} else {df$dotted<-df$max}
df$yrug<--(.2+(.1*(i-1)))
dflist[[i]]<-df
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
