taxdump=$1 #/path/to/taxdump.zip, must be in zip format
treefile=$2 #/path/to/ncbi.tre output
# e.g /path/to/database/megan/update-ncbi-tre.sh /path/to/database/taxonomy/taxdump/taxdump.zip /path/to/database/megan/ncbi.tre
/path/to/megan-ue/tools/ncbi/taxdmp2tree -i ${taxdump} -t ${treefile}
