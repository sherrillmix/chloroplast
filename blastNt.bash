#!/bin/bash
set -e

for ii in data/*trimmed.fastq;do
	outfile=work/$(basename $ii|sed 's/\.trimmed.fastq$//').blast.gz
	echo $ii -- $outfile
	awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' $ii|blastn -db ~/installs/db/nt/nt -num_threads 20 -max_target_seqs 20 -outfmt 6|gzip > $outfile
done
