#!/bin/bash
infile=$1
outfile=$2
[[ -z ${infile} ]] && exit

meganpath=/path/to/bin/megan/tools
dbpath=/path/to/database/megan
echo processing ${infile}
ln -sf ${dbpath}/ncbi.tre . && ln -sf ${dbpath}/ncbi.map .
paste <(${meganpath}/daa2info -i ${infile} -c2c Taxonomy | awk '{print $1}') <(${meganpath}/daa2info -i ${infile} -p -c2c Taxonomy | awk -F'\t' '{print $1,"\t",$2}') > ${outfile}_summary_count.tsv
