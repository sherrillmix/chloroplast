#!/bin/bash
set -e

for ii in data/*trimmed.fastq;do
	outfile=$(basename $ii|sed 's/\.trimmed.fastq$//').blast
	echo $ii -- $outfile
	blastn -query $ii  -db ~/installs/db/nt -num_threads 20 -max_target_seqs 20 -outfmt 6|gzip > $outfile
done
