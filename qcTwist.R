library(zoo)
library(LSD)
#install.packages("corrplot")
library(corrplot)

rootDir = "/home/labs/ulitsky/shared/data/tiling/Pum/run7/"
namesFile = "names.txt"
suffix = "combined.txt"
outFile = "qcTwist.pdf"
doCommonUMIs=FALSE

myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  rootDir = myArgs[[1]]
  namesFile = myArgs[[2]]
  suffix = myArgs[[3]]
  outFile = myArgs[[4]]
  if (length(myArgs)>4)
  {
    doCommonUMIs = as.logical(myArgs[[5]])
  }
}
names<-read.table(paste0(rootDir,namesFile),sep="\t",stringsAsFactors=FALSE,quote="",comment="",header=F)$V1

pdf(outFile,useDingbats = F,width=12)
tab<-read.table(paste0(rootDir,names[1],".",suffix),sep="\t",row.names=1)
tabs=list()
tabUMI<-data.frame(row.names = rownames(tab))
tabReads<-data.frame(row.names = rownames(tab))

for (name in sort(names))
{
  print(name)
  tab1<-read.table(paste0(rootDir,name,".",suffix),sep="\t",row.names=1,stringsAsFactors = F)
  tabUMI[[name]]=tab1$V2
  tabReads[[name]]=tab1$V3
  tabs[[name]]=tab1
}

par(mar=c(10,4,4,4))
barplot(apply(tabUMI,2,sum)/100000,las=2,ylab="Total UMIs (x100,000)",main="Total UMIs")
barplot(apply(tabReads,2,sum)/1000000,las=2,ylab="Total Reads (x1,000,000)",main="Total Reads")
barplot(apply(tabUMI,2,function(x){sum(x>0)/nrow(tabUMI)}),las=2,ylab="Represented tile %",main="Represented tiles",ylim=c(0,1))
abline(h=1)
barplot(apply(tabReads,2,function(x){sum(x>25)/nrow(tabReads)}),las=2,ylab="Represented tile %",main="Represented tiles with >=25 reads",ylim=c(0,1))
abline(h=1)
barplot(apply(tabUMI,2,function(x){sum(x>25)/nrow(tabUMI)}),las=2,ylab="Represented tile %",main="Represented tiles with >=25 UMIs",ylim=c(0,1))
abline(h=1)

boxplot(sapply(1:length(tabUMI),function(i){tabReads[[i]]/tabUMI[[i]]}),las=2,ylab="Reads/UMI",main="Reads/UMI",names=names(tabReads),pch=20)
abline(h=0)

tabGroups1 = table(unlist(lapply(strsplit(rownames(tabReads),"_"),"[",1)))
tabGroups2 = table(unlist(lapply(strsplit(rownames(tabReads),"_"),"[",2)))
tabGroups3 = table(unlist(lapply(strsplit(rownames(tabReads),":"),"[",1)))
tabGroups4 = table(unlist(lapply(strsplit(rownames(tabReads),":"),"[",2)))
tabGroups=c(tabGroups1,tabGroups2,tabGroups3,tabGroups4)
largeSubsets=names(tabGroups)[tabGroups>=20]

par(mar=c(4,4,4,4))
commonUMIshare=c()
p1UMIshare=c()
scaleFactors25 = c()
scaleFactors1 = c()
for (i in 1:length(tabUMI))
{
  par(mfrow=c(2,2))
  par(mar=c(4,4,4,4))
  nUMIs = tabUMI[[i]]
  print(paste("Plotting",names(tabUMI)[i]))
  heatscatter(tabReads[[i]],tabUMI[[i]],main=names(tabUMI)[i],cor=T,xlab="Reads",ylab="UMIs",method="spearman")
  abline(a=0,b=sum(tabUMI[[i]])/sum(tabReads[[i]]))
  heatscatter(tabUMI[[i]],tabReads[[i]]/tabUMI[[i]],xlab="UMIs",ylab="Reads/UMI",main="Reads/UMI")
  tMax=sapply(tabs[[i]]$V5,function(x){max(as.numeric(unlist(lapply(strsplit(unlist(strsplit(x,",")),"="),"[",2))))})
  heatscatter(tabUMI[tMax>0,i],tMax[tMax>0],main="Max(Reads)/UMI",xlab="UMI",ylab="Max reads/UMI")
  r1=rollapply((tabReads[[i]]/tabUMI[[i]])[order(tabUMI[[i]])],20,median,partial=T)
  plot(r1,pch=20,main="rolling median of read/UMI",xlab="UMIs",ylab="median reads/UMI in window")
  if (doCommonUMIs)
  {
    umiCounts = table(unlist(sapply(tabs[[i]]$V5,function(x){unlist(lapply(strsplit(unlist(strsplit(x,",")),"="),"[",1))})))
    plot(umiCounts[order(umiCounts)],pch=20)
    commonUMIshare[i]=sum(umiCounts[umiCounts>=100])/sum(umiCounts)
    p1UMIshare[i]=sum(umiCounts[umiCounts>=quantile(umiCounts,.99)])/sum(umiCounts)
  }
  myBreaks = seq(0,sum(nUMIs),10000)
  coverFun25 = sapply(myBreaks,function(x){sum((x/sum(nUMIs))*nUMIs>=25)})
  coverFun1 = sapply(myBreaks,function(x){sum((x/sum(nUMIs))*nUMIs>=1)})
  if (max(tabUMI[[i]])+1>500)
  {
	barplot(table(cut(tabUMI[[i]],include.lowest = T,breaks = c(-1,0,1,25,50,100,200,500,max(tabUMI[[i]])+1))),las=2,ylab="# tiles",main=paste(names(tabUMI)[i],"# of UMIs pre tile"))
  } else {
	barplot(table(cut(tabUMI[[i]],include.lowest = T,breaks = c(-1,0,1,25,50,100,200,500))),las=2,ylab="# tiles",main=paste(names(tabUMI)[i],"# of UMIs pre tile"))
  } 
  plot(myBreaks/10000,coverFun25,pch=20,las=2,xlab="nUMIs (x100,000)",ylab="Tiles with >=25 UMIs",main=paste(names(tabUMI)[i],"25X coverage saturation"))
  plot(myBreaks/10000,coverFun1,pch=20,las=2,xlab="nUMIs (x100,000)",ylab="Tiles with >=1 UMIs",main=paste(names(tabUMI)[i],"1X coverage saturation"))
  scaleFactors25=c(scaleFactors25,myBreaks[min(which(coverFun25>=max(coverFun25)*0.98))]/sum(nUMIs))
  scaleFactors1=c(scaleFactors1,myBreaks[min(which(coverFun1>=max(coverFun1)*0.98))]/sum(nUMIs))
  if (length(largeSubsets)>0)
  {
    par(mar=c(10,5,5,5))
    subCoverage= sapply(largeSubsets,function(x){sum(tabUMI[grep(x,rownames(tabUMI)),i]>=25)/length(grep(x,rownames(tabUMI)))})
    barplot(subCoverage,names=unlist(lapply(strsplit(names(subCoverage),"\\."),"[",1)),las=2,ylab="% of UMI coverage")
  }
}
meanReads = apply(tabReads,1,mean)
par(mfrow=c(1,1))
par(mar=c(12,5,5,5))
barplot(scaleFactors25*apply(tabReads,2,sum)/1000000,las=2,ylab="# reads to reach 98% of current 25X coverage (x1,000,000)",main=c("Coverage needed for 25X"))
barplot(scaleFactors1*apply(tabReads,2,sum)/1000000,las=2,ylab="# reads to reach 98% of current 1X coverage (x1,000,000)",main=c("Coverage needed for 1X"))
if (doCommonUMIs)
{
  par(mfrow=c(1,1))
  par(mar=c(10,4,4,4))
  barplot(100.0*commonUMIshare,names=names(tabUMI),ylab="% of UMIs that are seen in >100 tiles",main="CommonUMIs",las=2,ylim=c(0,100))
  barplot(100.0*p1UMIshare,names=names(tabUMI),ylab="% of UMIs that are in the 1% most common",main="CommonUMIs",las=2,ylim=c(0,100))
}
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
if (length(tabReads)>2)
{
  corrplot(cor(tabReads,method="spearman"),order="hclust",main="Raw read number")
  corrplot(cor(t(scale(t(tabReads[apply(tabReads,1,mean)>=10,]))),method="spearman"),order="hclust",main="Z scores of reads")
}
heatscatter(meanReads[meanReads>=10],apply(tabReads[meanReads>=10,],1,function(x){sd(x)/mean(x)}),xlab="Mean reads",ylab="CV reads",log="x",main="Data CV")

dev.off()
