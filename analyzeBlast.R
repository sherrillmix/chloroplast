
source('taxaId.R')

blastFiles<-list.files('work','.blast.gz$',full.names=TRUE)
taxas<-lapply(blastFiles,function(ii){
  message(ii)
  outFile<-sprintf('%s_taxa.csv',sub('.blast.gz$','',ii))
  if(file.exists(outFile)){
    message(' Reading ',outFile)
    taxaAssigns<-read.csv(outFile,row.names=1)
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
  taxa$best<-apply(taxa[,c('superkingdom','phylum','class','order','family','genus','species')],1,lastNotNa)
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
pdf('out/heat.pdf',width=10)
  heatmap(t(selectOtu),col=rev(heat.colors(100)),margins=c(10,4),scale='none')
dev.off()

animals<-sub('([A-Z]+[0-9]+).*','\\1',colnames(selectOtu))
selectOtu<-selectOtu[order(apply(selectOtu,1,sum)),order(animals)]
pdf('out/heatSort.pdf',width=9,height=13)
  #heatmap(selectOtu,col=rev(heat.colors(100)),margins=c(2,10),Colv=NA,scale='none')
  par(mar=c(.1,15,7,.1))
  image(1:ncol(selectOtu),1:nrow(selectOtu),t(selectOtu),col=rev(heat.colors(100)),xlab='',ylab='',xaxt='n',yaxt='n')
  axis(3,1:ncol(selectOtu),colnames(selectOtu),las=2)
  axis(2,1:nrow(selectOtu),rownames(selectOtu),las=1)
  abline(v=which(animals[-1]!=animals[-length(animals)])+.5)
  abline(h=2:nrow(selectOtu)-.5,col='#00000033')
  abline(v=2:ncol(selectOtu)-.5,col='#00000011')
  box()
dev.off()
