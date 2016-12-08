library(dnar)
library(levenR)
library(dnaplotr)
ALIGN<-FALSE

if(!exists('work/rbclAlign.blast.gz'))system('./alignRbcl.bash')

if(!exists('blast')){
  blast<-read.blast('work/rbclAlign.blast.gz')
  blast<-blast[blast$bit==ave(blast$bit,blast$qName,FUN=max),]
  blast<-blast[blast$score>500,]
  blast<-blast[ave(blast$qName,blast$qName,FUN=length)==1,]
  rownames(blast)<-blast$qName
}

ref<-read.fa('rbcl/refRbcl.fa')
if(!exists('seqs')){
  buffer<-400
  seqs<-read.fa('rbcl/rbcl.fa.gz')
  seqs$short<-sub(' .*','',seqs$name)
  seqs<-seqs[seqs$short %in% blast$qName,]
  seqs[,c('start','end')]<-blast[seqs$short,c('qStart','qEnd')]
  seqs$isPos<-blast[seqs$short,'tEnd']>blast[seqs$short,'tStart']
  #depending on substring not to care about negative start or ends past string
  seqs$trim<-substring(seqs$seq,seqs$start-buffer,seqs$end+buffer)
  seqs$trim[!seqs$isPos]<-revComp(seqs$trim[!seqs$isPos])
  write.fa(c('ref',seqs$short),c(ref$seq,seqs$trim),'work/trimRbcl.fa')
}
if(ALIGN){
  if(!file.exists('work/rbclAlign.fa')){
    system('mafft --retree 1 --maxiterate 0 --thread 50 work/trimRbcl.fa>work/rbclAlign.fa')
  }
  mafft<-read.fa('work/rbclAlign.fa')
  png('test.png',height=5000,width=5000,res=500)
    plotDNA(mafft$seq)
  dev.off()

  lAlign<-cacheOperation('work/rbclAlign.Rdat',levenAlign,seqs$trim,ref$seq,substring1=TRUE,nThread=56)
  dists<-leven(seqs$trim,ref$seq,substring1=TRUE,nThread=56)
  png('test.png',height=5000,width=5000,res=400)
    plotDNA(lAlign[[2]][order(dists)])
    abline(v=range(which(!strsplit(lAlign[[1]],'')[[1]] %in% c('.','-')))+c(-.5,.5))
  dev.off()
}



