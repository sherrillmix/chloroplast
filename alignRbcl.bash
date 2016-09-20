#!/bin/bash
set -e

makeblastdb -in rbcl/refRbcl.fa -dbtype nucl -out rbcl/rbclDb

zcat rbcl/rbcl.fa.gz |parallel --block-size 100M --recstart '>' -j 50 --pipe blastn -db rbcl/rbclDb -outfmt 6|gzip>work/rbclAlign.blast.gz

