#too much hassle. just download from ENA 
#http://www.ebi.ac.uk/ena/data/search?query=rbcl
library("rentrez")
library(parallel)
x<-entrez_search(db="nuccore", term="rbcl",retmax=as.integer(1e7),use_history=TRUE)
#seqs<-entrez_fetch('nuccore',web_history=x$web_history,rettype='fasta',retmax=as.integer(1e7))

stepSize<-1
allSeq<-mclapply(seq(1,length(x$ids),stepSize),function(ii){
  message(ii)
  end<-min(length(x$ids),ii+stepSize-1)
  seqs<-entrez_fetch('nuccore',web_history=x$web_history,rettype='fasta',retstart=ii,retmax=stepSize)
  if(length(gregexpr('^>|\n>',seqs)[[1]])!=end-ii+1)stop(simpleError('Missing sequences'))
  return(seqs)
},mc.cores=15)


