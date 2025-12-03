# MetaDetector

## Database Installation

### MEGAN
MEGAN is a tool for visualization of taxonomic outputs. This script downloads and installs the software and required databases.
```
bash download_MEGAN.sh
```

---
### NR and NT Databases
NR and NT are databases from NCBI that contain protein and nucleotide sequences, respectively.
```
bash download_NR_NT.sh
```

---
### Taxdump
The taxonomy file downloaded here are used to map taxonomic IDs in the NR and NT databases
```
bash download_taxdump.sh
```

---
### Contaminants
The databases downloaded in this script are used to remove common contaminants from your sequencing data.
```
bash download_contaminants.sh
```

---
#### Updating MEGAN
Updates databases used by MEGAN.
```
bash update_MEGAN_dbs.sh
```


