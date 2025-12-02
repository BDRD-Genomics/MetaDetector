#!/bin/bash
dbfile=${1}
protfile=${2}
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
EOF
