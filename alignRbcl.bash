#!/bin/bash
set -e

#https://www.ncbi.nlm.nih.gov/gene/844754
#https://www.ncbi.nlm.nih.gov/gene/844797
makeblastdb -in rbcl/refRbcl.fa -dbtype nucl -out rbcl/rbclDb
makeblastdb -in matk/refMatk.fa -dbtype nucl -out matk/matkDb

#http://www.ebi.ac.uk/ena/data/search?query=rbcl
#http://www.ebi.ac.uk/ena/data/search?query=matk
zcat rbcl/rbcl.fa.gz |parallel --block-size 100M --recstart '>' -j 50 --pipe blastn -db rbcl/rbclDb -outfmt 6|gzip>work/rbclAlign.blast.gz
zcat matk/matk.fa.gz |parallel --block-size 100M --recstart '>' -j 50 --pipe blastn -db matk/matkDb -outfmt 6|gzip>work/matkAlign.blast.gz
