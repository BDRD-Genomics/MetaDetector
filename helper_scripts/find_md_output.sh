#!/bin/bash
# e.g bash find_md_output.sh /full/path/to/search/dir file_name.tar.gz
finddir=$1
outfile=$2
#find ${finddir} \( -type f -wholename "*fastqc/posttrim/*_host_contaminant_rRNA_removed_pe_fastqc*" -o -wholename "*fastqc/pretrim/*_R1_fastqc.html" -o -wholename "*report*" -o -wholename "*pe_trim/contigs.fasta" -o -wholename "*pe_trim/scaffolds.fasta" -o -name "*_sorted.*" -o -wholename "*/blast/*.daa" \) | tar -cT - | pigz -9 -p 127 > ${outfile}
#find ${finddir} \( -type f -wholename "*/blast/*.daa" \) | tar -Ipigz -cT - | split -b 7G - ${outfile}
#find ${finddir} \( -type f -wholename "*/spades/meta_pe_trim/contigs.fasta" -o -wholename "*/blast/*" -o -wholename "*/mapping/_sorted.bam" \) | tar -Ipigz -cT - > ${outfile}
find ${finddir} \( -type f -wholename "*/spades/meta_pe_trim/contigs.fasta" -o -wholename "*/blast/*" \) | tar -Ipigz -cT - | split -b 2G - ${outfile}
