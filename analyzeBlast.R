
source('taxaId.R')

blastFiles<-list.files('work','.blast.gz$',full.names=TRUE)
taxas<-lapply(blastFiles[1:18],function(ii){
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

