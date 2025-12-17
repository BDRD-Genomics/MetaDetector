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
scripts=/path/to/update_scripts

bash $scripts/update_megan-map-db.sh $md_db_dir/megan/megan-map-$(date +%d%b%y).mdb prot.accession2taxid
bash $scripts/update_megan-nucl-db.sh $md_db_dir/megan/megan-nucl-$(date +%d%b%y).mdb prot.accession2taxid



