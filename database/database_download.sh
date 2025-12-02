#############################################################################
# Author: Gregory Rice <gregory.k.rice.ctr@health.mil>
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

# Modify paths
md_db_dir=/path/to/databases
scripts=/path/to/helper_scripts

# Requirements include: nr, nt and megan

# Optional automated download: contaminants, rRNA and human
cd $md_db_dir/ref_genomes
mkdir human silva
wget https://edge-dl.lanl.gov/EDGE/light/human_ref_GRCh38_all.fa.gz
wget -cvb -O $md_db_dir/ https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz
wget -cvb -O $md_db_dir/ https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
zcat SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz > SILVA_138.2_LSU-SSURef_NR99_tax_silva.fasta && pigz SILVA_138.2_LSU-SSURef_NR99_tax_silva.fasta
wget -cvb -O $md_db_dir/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz

# Download and install MEGAN7
# only download this once
mkdir $md_db_dir/megan
wget -cvb -O $md_db_dir/megan/MEGAN_Community_unix_7_1_1.sh https://software-ab.cs.uni-tuebingen.de/download/megan7/MEGAN_Community_unix_7_1_1.sh
bash $md_db_dir/megan/MEGAN_Community_unix_7_1_1.sh  # Respond to prompts for installation
wget -cvb -O $md_db_dir/megan/megan-nr-r2.zip https://software-ab.cs.uni-tuebingen.de/download/megan7/megan-nr-r2.zip
unzip -d $md_db_dir/megan/ $md_db_dir/megan/megan-nr-r2.zip && mv $md_db_dir/megan/megan-nr-r2.mdb $md_db_dir/megan/megan-map-$(date +%d%b%y).mdb
ln -sf $md_db_dir/megan/megan-map-$(date +%d%b%y).mdb $md_db_dir/megan/megan-map.mdb
wget -cvb -O $md_db_dir/megan/megan-genome-r1.zip https://software-ab.cs.uni-tuebingen.de/download/megan7/megan-genome-r1.zip
unzip -d $md_db_dir/megan/ $md_db_dir/megan/megan-genome-r1.zip && mv $md_db_dir/megan/megan-genome-r1.mdb $md_db_dir/megan/megan-nucl-$(date +%d%b%y).mdb
ln -sf $md_db_dir/megan/megan-nucl-$(date +%d%b%y).mdb $md_db_dir/megan/megan-nucl.mdb


# Update NR and NT Commands:
update_blastdb.pl --decompress nr && mv nr $md_db_dir/nr_$(date +%d%b%y)
ln -sf $md_db_dir/nr_$(date +%d%b%y) $md_db_dir/nr
update_blastdb.pl --decompress nt && mv nt $md_db_dir/nt_$(date +%d%b%y)
ln -sf $md_db_dir/nt_$(date +%d%b%y) $md_db_dir/nt
 
# Creating Diamond Databases: 
mkdir $md_db_dir/nr_dmnd_$(date +%d%b%y) && ln -sf $md_db_dir/nr_dmnd_$(date +%d%b%y) $md_db_dir/nr && cd $md_db_dir/nr
blastdbcmd blastdbcmd -db $md_db_dir/nr/nr -entry all | pigz -9 > $md_db_dir/nr_dmnd/nr.gz
cd $md_db_dir/nr_dmnd
diamond makedb --in nr.gz --db nr.dmnd

# Create vhunter_acc.db:
# Redownload this to update every once in awhile
mkdir $md_db_dir/taxdump_$(date +%d%b%y) && ln -sf $md_db_dir/taxdump_$(date +%d%b%y) $md_db_dir/taxdump
wget -cvb -O $md_db_dir/taxdump/nucl_gb.accession2taxid.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz #(or most recent nucl taxid)
wget -cvb -O $md_db_dir/taxdump/prot.accession2taxid.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz #(or most recent prot taxid)
cd $md_db_dir/taxdump/
pigz -dc nucl_gb.accession2taxid.gz > nucl_gb.accession2taxid
pigz -dc nucl_gb.accession2taxid.gz > prot.accession2taxid

bash $scripts/update_megan-map.db $md_db_dir/megan/megan-map-$(date +%d%b%y).mdb prot.accession2taxid
bash $scripts/update_megan-nucl.db $md_db_dir/megan/megan-nucl-$(date +%d%b%y).mdb prot.accession2taxid


