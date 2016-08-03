library("rentrez")
x<-entrez_search(db="nuccore", term="rbcl",retmax=as.integer(1e7))
#seqs<-entrez_fetch('nuccore',x$ids,rettype='fasta')

allSeq<-c()
stepSize<-10
for(ii in seq(1,length(x$ids),stepSize)){
  message(ii)
  end<-max(length(x$ids),ii+stepSize-1)
  seqs<-entrez_fetch('nuccore',x$ids[ii:end],rettype='fasta')
  if(length(gregexpr('>',seqs)[[1]])!=end-ii+1)stop(simpleError('Missing sequences'))
  allSeq<-c(allSeq,seqs)
  Sys.sleep(2)
}

#while(length(tmp<-entrez_search(db="nuccore", term="rbcl",retstart=as.integer(pos),retmax=as.integer(blockSize))$ids)>0){
  #message(pos)
  #ids<-c(ids,tmp)
  #pos<-pos+blockSize
  #message('Sleep')
  #Sys.sleep(1)
  #if(potentialProblem)stop(simpleError('Incomplete pull not final pull'))
  #if(length(tmp)!=blockSize)potentialProblem<-TRUE
#}

