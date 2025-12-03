# MetaDetector

## Database Installation

Database installation overview. Multiple scripts, here is what they do.

### MEGAN
MEGAN is a tool for visualization of taxonomic analysis. This script downloads and installs the software and required databases.

bash download_MEGAN.sh

### NR and NT Databases
NR and NT are databases from NCBI that contain protein and nucleotide sequences, respectively.

bash download_NR_NT.sh


### DIAMOND Databases
Diamond is a tool for accelerated BLAST jobs. This step that the NR database is listed above is fully downloaded.

bash create_DIAMOND_db.sh


### Download Taxdump
The taxonomy file downloaded here are used to map taxonomic IDs in the NR and NT databases

bash download_taxdump.sh

### Contaminants
The databases downloaded in this script are used to remove common contaminants from your sequencing data.

bash download_contaminants.sh


#### Updating MEGAN
Updates databases used by MEGAN.

bash update_MEGAN_dbs.sh



