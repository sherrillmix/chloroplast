library("rentrez")
x<-entrez_search(db="nuccore", term="rbcl",retmax=as.integer(1e7))
#seqs<-entrez_fetch('nuccore',x$ids,rettype='fasta')

#don't need this. just need to make sure retmax is an integer
#pos<-0
#ids<-c()
#potentialProblem<-FALSE
#blockSize<-10000
#while(length(tmp<-entrez_search(db="nuccore", term="rbcl",retstart=as.integer(pos),retmax=as.integer(blockSize))$ids)>0){
  #message(pos)
  #ids<-c(ids,tmp)
  #pos<-pos+blockSize
  #message('Sleep')
  #Sys.sleep(1)
  #if(potentialProblem)stop(simpleError('Incomplete pull not final pull'))
  #if(length(tmp)!=blockSize)potentialProblem<-TRUE
#}

