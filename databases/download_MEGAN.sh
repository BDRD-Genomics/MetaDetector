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


