
library(dnar)
source('taxaId.R')
targetTaxa<-c('superkingdom','phylum','class','order','family','genus','species')

blastFiles<-list.files('work','.blast.gz$',full.names=TRUE)
taxas<-lapply(blastFiles,function(ii){
  message(ii)
  outFile<-sprintf('%s_taxa.csv',sub('.blast.gz$','',ii))
  if(file.exists(outFile)){
    message(' Reading ',outFile)
    taxaAssigns<-read.csv(outFile,row.names=1,stringsAsFactors=FALSE)
  }else{
    message(' Creating ',outFile)
    x<-read.table(ii,stringsAsFactors=FALSE)
    colnames(x)<-c('qName','tName','percIdent','alignLength','mismatch','gapOpens','qStart','qEnd','tStart','tEnd','eVal','bitScore')
    x$accession<-sapply(strsplit(x$tName,'\\|'),'[[',4)
    x$taxa<-accessionToTaxa(x$accession,'dump/accessionTaxa.sql')
    x$maxBit<-ave(x$bitScore,x$qName,FUN=max)
    x<-x[x$bitScore==x$maxBit&!is.na(x$taxa),]
    taxonomy<-getTaxonomy(x$taxa,taxaNodes,taxaNames)
    taxaAssigns<-do.call(rbind,by(as.data.frame(taxonomy,stringsAsFactors=FALSE),x$qName,FUN=condenseTaxa))
    bestScore<-x[!duplicated(x$qName),c('qName','alignLength','percIdent')]
    rownames(bestScore)<-bestScore$qName
    taxaAssigns<-cbind(taxaAssigns,bestScore[rownames(taxaAssigns),c('alignLength','percIdent')])
    write.csv(taxaAssigns,outFile)
  }
  return(taxaAssigns)
})
names(taxas)<-basename(blastFiles)
taxas<-mclapply(taxas,function(taxa){
  cat('.')
  taxa$best<-apply(taxa[,targetTaxa],1,lastNotNa)
  return(taxa)
},mc.cores=10)

taxaTab<-lapply(taxas,function(x)sort(table(x$best)))
topTaxa<-sapply(taxaTab,function(x)tail(x/sum(x),1))
propPlasmodium<-sapply(taxas,function(x)mean(x$order %in% 'Haemosporida'))
topDf<-data.frame('file'=names(taxas),'topTaxa'=names(topTaxa),'topProportion'=round(topTaxa,3),'haemosporidaProportion'=round(propPlasmodium,3),stringsAsFactors=FALSE)
write.csv(topDf,'out/topSpecies.csv',row.names=FALSE)

isMalaria<-grepl('cytb1',names(taxaTab))
taxaTab<-taxaTab[!isMalaria]
taxas<-taxas[!isMalaria]
#pooling R1 and R2
otuTab<-table(unlist(lapply(taxas,function(x)x$best)),rep(sub('_R[12]_','_',names(taxas)),sapply(taxas,nrow)))

selectOtu<-apply(otuTab,2,function(x)x/sum(x))
selectOtu<-selectOtu[apply(selectOtu,1,max)>.01,]
colnames(selectOtu)<-sub('_S[0-9]+_L001_001.blast.gz','',colnames(selectOtu))
animals<-sub('([A-Z]+[0-9]+).*','\\1',colnames(selectOtu))
selectOtu<-selectOtu[order(apply(selectOtu,1,sum)),order(animals)]

insetLegend<-function(breaks,cols,insetPos=c(.015,.025,.25,.04)){
  insetPos<-c(grconvertX(insetPos[1],'nfc','user'),grconvertY(insetPos[2],'nfc','user'),grconvertX(insetPos[3],'nfc','user'),grconvertY(insetPos[4],'nfc','user'))
  breakPos<-((breaks[-1])-(min(breaks[-1])))/max((breaks[-1])-(min(breaks[-1])))*(insetPos[3]-insetPos[1])+insetPos[1]
  rect(breakPos[-1]+1e-3,insetPos[2],breakPos[-length(breakPos)],insetPos[4],col=cols[-1],xpd=NA,border=NA)
  rect(insetPos[1],insetPos[2],insetPos[3],insetPos[4],xpd=NA)
  prettyLabs<-pretty(breaks)
  prettyLabs<-prettyLabs[prettyLabs<max(breaks)]
  prettyPos<-prettyLabs
  prettyPos<-(prettyLabs-(min(breaks[-1])))/((max(breaks[-1]))-(min(breaks[-1])))*(insetPos[3]-insetPos[1])+insetPos[1]
  segments(prettyPos,insetPos[2],prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.1,xpd=NA)
  text(prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.175,prettyLabs,xpd=NA,adj=c(.5,1),cex=.85)
  text(mean(insetPos[c(1,3)]),insetPos[4]+diff(insetPos[c(2,4)])*.45,"Proportion of sample",xpd=NA,adj=c(.5,0))
}

breaks<-c(0,seq(min(selectOtu),max(selectOtu),length.out=501))
cols<-c('white',tail(rev(heat.colors(520)),500))
pdf('out/heatCluster.pdf',width=10)
  heatmap(t(selectOtu),col=cols,breaks=breaks,margins=c(10,4),scale='none',Colv=NA)
  insetLegend(breaks,cols)
dev.off()

pdf('out/heatSort.pdf',width=9,height=13)
  #heatmap(selectOtu,col=rev(heat.colors(100)),margins=c(2,10),Colv=NA,scale='none')
  par(mar=c(.1,15,7,.1))
  image(1:ncol(selectOtu),1:nrow(selectOtu),t(selectOtu),col=cols,xlab='',ylab='',xaxt='n',yaxt='n')
  axis(3,1:ncol(selectOtu),colnames(selectOtu),las=2)
  axis(2,1:nrow(selectOtu),rownames(selectOtu),las=1)
  abline(v=which(animals[-1]!=animals[-length(animals)])+.5)
  abline(h=2:nrow(selectOtu)-.5,col='#00000033')
  abline(v=2:ncol(selectOtu)-.5,col='#00000011')
  box()
  insetLegend(breaks,cols,insetPos=c(.05,.93,.25,.945))
dev.off()

sumId<-do.call(rbind,mclapply(taxas,function(xx){
  #need to fill left on this for skipped taxa steps
  xx<-apply(xx[,targetTaxa],1,function(yy){
    rev(fillDown(rev(yy),errorIfFirstEmpty=FALSE))
  })
  superkingdom<-sum(!is.na(xx[1,]))
  xx[,xx[1,]!='Eukaryota']<-NA
  onlyEuk<-apply(xx,1,function(yy)sum(!is.na(yy)))
  return(c(superkingdom,onlyEuk))
},mc.cores=10))
fastqFiles<-file.path('data',sub('.blast.gz','.trimmed.fastq',names(taxas)))
readCounts<-as.numeric(cacheOperation('work/readCounts.Rdat',mclapply,fastqFiles,function(x)system(sprintf("wc -l %s|cut -f1 -d' '",x),intern=TRUE),mc.cores=10))/4

propId<-sumId/matrix(readCounts,nrow=nrow(sumId),ncol=ncol(sumId))
colnames(propId)<-c(targetTaxa[1],'Eukaryota',targetTaxa[-1])

primers<-sub('KSG[0-9]+','',sub('_.*','',rownames(propId)))
cols<-rainbow.lab(length(unique(primers)),alpha=.6)
names(cols)<-unique(primers)

pdf('out/mapped.pdf')
par(mfrow=c(2,2),mar=c(4,4,1,.1))
for(primer in unique(primers)){
  plot(1,1,type='n',ylab='Proportion of reads mapped',xlab='',xaxt='n',ylim=c(0,1),las=1,xlim=c(1,ncol(propId)),main=primer,mgp=c(2.5,1,0))
  #axis(1,1:ncol(propId),sub('super','',colnames(propId)))
  axis(1,1:ncol(propId),FALSE)
  text(1:ncol(propId), y=convertLineToUser(.8), labels=sub('super','',colnames(propId)), srt=45, adj=c(1,1), xpd=TRUE)
  thisProps<-propId[primers==primer,]
  paired<-aggregate(thisProps,list(sub('R[12]_.*','',rownames(thisProps))),function(x){
    if(abs(log2(x[1]/x[2]))>log2(1.2))stop(simpleError('Mismatched pair'))
    mean(x)
  })[,-1] #remove group column
  apply(paired,1,function(x)lines(1:ncol(propId),x,col='#00000077',lwd=2))
}
dev.off()

kingTab<-table(unlist(lapply(taxas,function(x)x$superkingdom)),rep(names(taxas),sapply(taxas,nrow)))


names(taxas)
baseName<-sub('_S[0-9]+_L001_R[12]_.*$','',names(taxas))

rares<-lapply(unique(baseName),function(xx){
  thisTaxas<-unlist(lapply(taxas[baseName==xx],function(yy)yy$best))
  return(rareEquation(table(thisTaxas),round(seq(.01,1,.005)*length(thisTaxas))))
})
names(rares)<-unique(baseName)
primers<-sub('KSG[0-9]+','',sub('_.*','',names(rares)))
xlim<-c(0,max(as.numeric(unlist(lapply(rares,names)))))
ylim<-c(1,max(as.numeric(unlist(rares))))
pdf('out/rarefaction.pdf')
par(mfrow=c(2,2),mar=c(4,4.3,1,.3))
for(primer in unique(primers)){
  thisRares<-rares[primers==primer]
  plot(1,1,type='n',ylab='Number of taxa',xlab='',ylim=ylim,las=1,xlim=xlim/100000,main=primer,mgp=c(3,1,0))
  title(xlab='Reads sequenced (100,000)',mgp=c(2,1,0))
  sapply(thisRares,function(x)lines(as.numeric(names(x))/100000,x,col='#00000077',lwd=2))
}
dev.off()


