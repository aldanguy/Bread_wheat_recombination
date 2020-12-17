#!/bin/bash

## Define scaffolds limits of v1.0 assembly

base=${1}
source ${base}




echo -e "chr\tposlscaf\tposrscaf" > ${r_scaffolds}scaffolds.txt

<<COMMENTS
cat ${scaffolds_input} |
grep -e"^chr[0-9]" |
awk '$6 != "gap" { print $0 }'|
cut -f1,2,3 |
sed "s/chr//" |
sort -b -V -k1,1 -k2,2 |
sed 's| |\t|g'>> ${r_sources}scaffolds.txt

COMMENTS

cat ${scaffolds_input} |
grep -e"^chr[0-9]" |
awk '$6 != "gap" { print $0 }'|
cut -f1,2,3 |
sed "s/chr//" |
sort -b -k1,1 -k2,2 |
awk 'BEGIN {FS="\t"} NR==1 {chr=$1; poslscaf=$2; posrscaf=$3; next} ($1==chr && $2==posrscaf+1) {posrscaf=$3; next} {print chr"\t"poslscaf"\t"posrscaf; chr=$1; poslscaf=$2; posrscaf=$3} END {print chr"\t"poslscaf"\t"posrscaf}'|
sed 's| |\t|g' >> ${r_scaffolds}scaffolds.txt
