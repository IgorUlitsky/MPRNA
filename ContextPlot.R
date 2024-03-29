
myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  configFile = myArgs[[1]]
}

config = read.table(configFile,row.names = 1,stringsAsFactors = F,fill=T)

rootDir = config["rootDir",1]
File = config["finalTab",1]
condition  = config["condition",1]
numOfReplicates = config["numOfReplicates",1]
outFile = config["outFile",1]

#rootDir = "/home/labs/ulitsky/igoru/"
#File = "twistGFP_rep3.processTwist.finalTab.txt"
#condition="Cprime_"
#numOfReplicates=3

finalTab <- read.table(paste0(rootDir,File),sep="\t",stringsAsFactors=FALSE,quote="",comment="",header=T)

# Calculate median of replicates 
finalTab[,paste0(condition,"_med")] <- apply(finalTab[grep(condition, colnames(finalTab))], 1, median, na.rm = T)

## function
makePlot<-function(gene,subMat,windowSize,ylab,drawLegend=T,title=gene){
  ncA = subMat[,paste0(condition,"_med")]
  
  # calculate minimum and maximum values (for error bars)
  segA = subMat[,grep(condition, colnames(subMat))][1:numOfReplicates]
  segA$min <- apply(segA, 1, min, na.rm=T)
  segA$max <- apply(segA, 1, max, na.rm=T)
  segA$Index <- subMat[,"Index"]
  #segA <- segA[order(segA$Index),]
  
  minAll = min(c(ncA),na.rm=T)
  maxAll = max(c(ncA),na.rm=T)
  
  usable=rep(TRUE,length(ncA))
  myOrder <- order(subMat$Index)[usable]
  #myCols = sapply(subMat$Repetitive,function(x){ if (x) {"black"} else {"darkgrey"}})[myOrder]
  
  for (myStart in seq(1,length(ncA),windowSize)){
    myEnd = min(myStart+windowSize,length(ncA))
    if (!is.na(minAll)){
      myTitle = title
      if (myStart>1 || myEnd <length(ncA)){
        myTitle  = paste0(title," ",myStart,":",(myEnd-1))
      }
      plot(myStart:myEnd,ncA[myOrder[myStart:myEnd]],main=myTitle,ylim=c(min(-1.5,minAll),max(1.5,maxAll)),col="red",pch=20,ylab=ylab,xlab="Tile number",bty="n")
      lines(myStart:myEnd,ncA[myOrder[myStart:myEnd]],type="l",col="red")
      for (i in myStart:myEnd){
        segments(x0 = i, 
                 y0 = segA[myOrder[i],"min"], 
                 y1 = segA[myOrder[i],"max"], col = "lightgray")
      }
      
      abline(h=0, col = "lightgray")
      
      if (drawLegend){
        legend("topright",c(condition),col=c("red"),pch=c(20),bg = "gray94",ncol = 2, x.intersp = 0.5, horiz = T)
      }
    }
  }
}


## generating context plots
windowSize= 100
par(mar=c(4,4,4,4))
par(mfrow=c(1,1))
pdf(outFile,width=12)
for (gene in unique(finalTab$GeneName)){
  print(paste(gene))
  subMat<-finalTab[grepl(paste0("^",gene,"$"),finalTab$GeneName),]
  makePlot(gene,subMat,windowSize,"Nuc/Cyto (log2)")
}
dev.off()

