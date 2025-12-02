#!/bin/bash
dbfile=${1} # /path/to/megan-map.db
protfile=${2} # /path/to/prot.accession2taxid
sqlite3 -batch ${dbfile} <<EOF
DROP TABLE IF EXISTS acc2taxid;
CREATE TABLE acc2taxid (
  accession TEXT PRIMARY KEY,
  accession_version UNIQUE,
  taxid integer NOT NULL,
  gi integer
);
.separator "\t"
.import ${protfile} acc2taxid
DELETE FROM acc2taxid WHERE accession_version="accession.version";
UPDATE mappings
SET Taxonomy = (SELECT taxid from acc2taxid where acc2taxid.accession = mappings.Accession) WHERE EXISTS (SELECT * from acc2taxid where acc2taxid.accession = mappings.Accession);
INSERT OR IGNORE INTO mappings (Accession,Taxonomy)
SELECT accession, taxid
FROM acc2taxid;
DROP TABLE IF EXISTS acc2taxid;
VACUUM;
EOF
