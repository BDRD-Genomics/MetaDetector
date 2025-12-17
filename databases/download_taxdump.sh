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

# Redownload this to update every once in awhile
mkdir $md_db_dir/taxdump_$(date +%d%b%y) && ln -sf $md_db_dir/taxdump_$(date +%d%b%y) $md_db_dir/taxdump
wget -cvb -O $md_db_dir/taxdump/nucl_gb.accession2taxid.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz #(or most recent nucl taxid)
wget -cvb -O $md_db_dir/taxdump/prot.accession2taxid.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz #(or most recent prot taxid)
cd $md_db_dir/taxdump/
pigz -dc nucl_gb.accession2taxid.gz > nucl_gb.accession2taxid
pigz -dc nucl_gb.accession2taxid.gz > prot.accession2taxid


