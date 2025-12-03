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

md_db_dir=/path/to/databases

# Creating Diamond Databases: 
mkdir $md_db_dir/nr_dmnd_$(date +%d%b%y) && ln -sf $md_db_dir/nr_dmnd_$(date +%d%b%y) $md_db_dir/nr && cd $md_db_dir/nr
blastdbcmd blastdbcmd -db $md_db_dir/nr/nr -entry all | pigz -9 > $md_db_dir/nr_dmnd/nr.gz
cd $md_db_dir/nr_dmnd
diamond makedb --in nr.gz --db nr.dmnd


