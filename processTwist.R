library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(zoo)
library(LSD)
library(corrplot)

configFile = "/home/labs/ulitsky/shared/data/mouse/CytoNuc/MN/Nzip2/config.txt"

myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  configFile = myArgs[[1]]
}

config = read.table(configFile,row.names = 1,stringsAsFactors = F,fill=T)

repeatFile=NA
featuresFile=NA

useUMIs = F
useReps = T
useDESeq2 = T

rootDir = config["rootDir",1]
groupsFile = config["groupsFile",1]
ratiosFile = config["ratiosFile",1]
groupRatiosFile = config["groupRatiosFile",1]
suffix = config["suffix",1]
outBase = config["outBase",1]
if ("repeatFile" %in% rownames(config)) repeatFile = config["repeatFile",1]
if ("featuresFile" %in% rownames(config)) featuresFile = config["featuresFile",1]
minWCEReads = as.numeric(config["minWCEReads",1])
minInputReads = as.numeric(config["minInputReads",1])
minRatioReads=as.numeric(config["minRatioReads",1])
pseudo=as.numeric(config["pseudo",1])
if ("useUMIs" %in% rownames(config)) useUMIs=as.logical(config["useUMIs",1])
if ("useReps" %in% rownames(config)) useReps=as.logical(config["useReps",1])
if ("useDESeq2" %in% rownames(config)) useDESeq2=as.logical(config["useDESeq2",1])

groups<-read.table(paste0(rootDir,groupsFile),sep="\t",stringsAsFactors=FALSE,quote="",comment="",header=F,fill=T)
names = groups$V1
ratios<-read.table(paste0(rootDir,ratiosFile),sep="\t",stringsAsFactors=FALSE,quote="",comment="",header=F)
groupRatios<-c()
if (file.size(paste0(rootDir,groupRatiosFile))>0)
  groupRatios<-read.table(paste0(rootDir,groupRatiosFile),sep="\t",stringsAsFactors=FALSE,quote="",comment="",header=F)

pdf(paste0(outBase,".pdf"),useDingbats = F,width=12)
tab<-read.table(paste0(rootDir,names[1],".",suffix),sep="\t",row.names=1)
tabs=list()
tabUMI<-data.frame(row.names = rownames(tab))
tabReads<-data.frame(row.names = rownames(tab))


for (name in names)
{
  print(name)
  tab1<-read.table(paste0(rootDir,name,".",suffix),sep="\t",row.names=1,stringsAsFactors = F)
  tabUMI[[name]]=tab1$V2
  tabReads[[name]]=tab1$V3
  tabs[[name]]=tab1
}

if (useUMIs)
{
  tab=tabUMI
} else {
  tab=tabReads
}
  

filteredTab = tab
print(paste("Pre-filtering:",nrow(tab)," tiles"))
if (length(grep("WCE",groups$V2))>0)
{
  wceAv = apply(cbind(tab[,grep("WCE",groups$V2)]),1,mean)
  filteredTab <- subset(tab,wceAv>=minWCEReads)
  print(paste("WCE filter:",nrow(filteredTab),"tiles",sprintf("%3.2f",100.0*nrow(filteredTab)/nrow(tab)),"%"))
}
if (length(grep("Input",groups$V2))>0)
{
  inputAv = apply(cbind(filteredTab[,grep("Input",groups$V2)]),1,mean)
  filteredTab <- subset(filteredTab,inputAv>=minInputReads)
  print(paste("Input filter:",nrow(filteredTab),"tiles",sprintf("%3.2f",100.0*nrow(filteredTab)/nrow(tab)),"%"))
}

norm = 1000000/apply(filteredTab,2,sum)
normTab<-data.frame(t(apply(filteredTab,1,function(x){x*norm})))

names(normTab)<-names(filteredTab)
row.names(normTab)<-row.names(filteredTab)


ratioTab = data.frame(sapply(1:nrow(ratios),function(i){print(i);x=log((pseudo+normTab[[ratios[i,2]]])/(pseudo+normTab[[ratios[i,3]]]),base=2);print(length(x));x}))
names(ratioTab)=ratios[,1]
for (i in 1:length(ratioTab))
{
  ratioTab[filteredTab[[ratios[i,2]]]+filteredTab[[ratios[i,3]]]<minRatioReads*2,i]
}

row.names(ratioTab)<-row.names(filteredTab)


ratioTab[["Index"]]<-as.numeric(unlist(lapply(strsplit(row.names(filteredTab),"_"),function(x){x[[2]]})))
#ratioTab[["SegmentName"]]<-unlist(lapply(strsplit(row.names(ratioTab),"_"),function(x){paste(x[1],x[2],sep="_")}))
ratioTab[["SegmentName"]]<-row.names(ratioTab)
ratioTab[["SegmentName"]][grep("BORG|FIRRE",ratioTab$SegmentName)]<-unlist(lapply(strsplit(row.names(ratioTab)[grep("BORG|FIRRE",ratioTab$SegmentName)],"_"),function(x){paste(x[1:3],collapse="_")}))
ratioTab[["SegmentName"]][grep("Intron",ratioTab$SegmentName)]<-unlist(lapply(strsplit(row.names(ratioTab)[grep("Intron",ratioTab$SegmentName)],"_"),function(x){paste(x[1:5],collapse="_")}))
ratioTab[["SegmentName"]][grep("Hybrid",ratioTab$SegmentName)]<-unlist(lapply(strsplit(row.names(ratioTab)[grep("Hybrid",ratioTab$SegmentName)],"_"),function(x){paste(x[1:4],collapse="_")}))
ratioTab[["GeneName"]]<-unlist(lapply(strsplit(row.names(ratioTab),"_"),function(x){x[1]}))
ratioTab[["WT"]][grep(":",rownames(ratioTab),invert=TRUE)]=T
ratioTab[["GeneName"]]<-unlist(lapply(strsplit(row.names(ratioTab),"_"),function(x){x[1]}))
if (!is.na(repeatFile))
{
  repeatInfo<-read.table(repeatFile,sep="\t",stringsAsFactors=FALSE,quote="",comment="")
  names(repeatInfo)<-c("Segment","Repeat","Score")
  ratioTab[["Repetitive"]]<-ratioTab$SegmentName %in% unique(repeatInfo$Segment)
  mergedTab=merge(ratioTab,aggregate(repeatInfo$Repeat,by=list(repeatInfo$Segment),FUN=paste),by.x="SegmentName",by.y="Group.1",all.x=TRUE)
} else {
  mergedTab=ratioTab
}

if (!is.na(featuresFile))
{
  segmentFeatures<-read.table(featuresFile,sep="\t",stringsAsFactors=FALSE,quote="",comment="",header=TRUE)
  segmentFeatures<-segmentFeatures[!grepl("_S1$",segmentFeatures$Name),]
  segmentFeatures$Name=unlist(lapply(strsplit(segmentFeatures$Name,"_"),function(x){paste(x[1],x[2],sep="_")}))
  mergedTab=merge(mergedTab,segmentFeatures,by.x="SegmentName",by.y="Name",all.x=TRUE)
}

names(filteredTab)=paste0(names(filteredTab),".raw")
mergedTab=merge(mergedTab,filteredTab,by.x="SegmentName",by.y=0,all.x=TRUE)
if (!is.na(repeatFile))
{
  mergedTab[["CombinedRepeats"]]<-as.character(unlist(lapply(mergedTab[,"x"],FUN=function(y){toString(unlist(y))})))
}

charCol<-unlist(mergedTab$x)
if (length(grep("^x$",names(mergedTab)))>0)
{
  finalTab = mergedTab[,-grep("^x$",names(mergedTab))]
} else 
{
  finalTab = mergedTab 
}
rownames(finalTab)=finalTab$SegmentName

if (useDESeq2 & sum(as.logical(groups$V4))>0)
{
  coldata = data.frame(condition=factor(groups[groups$V4,2]),replicate=factor(groups[groups$V4,3]))
  if (useReps)
  {
    ddsMat <- DESeqDataSetFromMatrix(countData = filteredTab[,groups$V4],
                                     colData = coldata,
                                     design = ~ condition+replicate)
  } else {
    ddsMat <- DESeqDataSetFromMatrix(countData = filteredTab[,groups$V4],
                                     colData = coldata,
                                     design = ~ condition)
  }
  
  dds <- DESeq(ddsMat)
  
  for (i in 1:nrow(groupRatios))
  {
    resC <- results(dds, alpha=0.05,contrast=c("condition",groupRatios[i,2],groupRatios[i,3]))
    plot(resC$log2FoldChange,-log(resC$pvalue,10),pch=20,col=1+(resC$padj<0.05),main=groupRatios[i,1])
    finalTab[[paste0(groupRatios[i,1],".log2FoldChange")]]=resC[finalTab$SegmentName,"log2FoldChange"]
    finalTab[[paste0(groupRatios[i,1],".pvalue")]]=resC[finalTab$SegmentName,"pvalue"]
    finalTab[[paste0(groupRatios[i,1],".padj")]]=resC[finalTab$SegmentName,"padj"]
    finalTab[[paste0(groupRatios[i,1],".baseMean")]]=resC[finalTab$SegmentName,"baseMean"]
  }
  if (length(grep("log2FoldChange",names(finalTab)))>1)
  {
    heatpairs(data.matrix(finalTab[,grep("log2FoldChange",names(finalTab))]))
  }
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  # Clustering of conditions
  rld <- vst(dds, blind=TRUE)
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$condition, rld$replicate,sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  par(mar=c(10,10,20,20))
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors,width =4,height =4,show_colnames = T,show_rownames = T)
  plotPCA(rld, intgroup=c("condition"))
}

write.table(finalTab,paste0(outBase,".finalTab.txt"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE,na="")

dev.off()
