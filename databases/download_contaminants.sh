#############################################################################
# Contact: Reachback Support <usn.detrick.nmrc.mbx.genomics-reach-back@health.mil>
# Additional Contact: Regina Cer <regina.z.cer.civ@health.mil>
#
# License:
# Metadetector is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, version X of the License or any later version. For
# details, please refer to https://www.gnu.org/licenses/
#############################################################################

eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)
conda activate md

md_db_dir=/path/to/databases

cd $md_db_dir
mkdir human silva

wget -cvb -O $md_db_dir/human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

wget -cvb -O $md_db_dir/silva/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz

wget -cvb -O $md_db_dir/silva/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz

zcat $md_db_dir/silva/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz $md_db_dir/silva/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz > $md_db_dir/silva/SILVA_138.2_LSU-SSURef_NR99_tax_silva.fasta && pigz $md_db_dir/silva/SILVA_138.2_LSU-SSURef_NR99_tax_silva.fasta


