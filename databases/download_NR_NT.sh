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

# Update NR and NT Commands:
update_blastdb.pl --decompress nr && mv nr $md_db_dir/nr_$(date +%d%b%y)
ln -sf $md_db_dir/nr_$(date +%d%b%y) $md_db_dir/nr
update_blastdb.pl --decompress nt && mv nt $md_db_dir/nt_$(date +%d%b%y)
ln -sf $md_db_dir/nt_$(date +%d%b%y) $md_db_dir/nt


