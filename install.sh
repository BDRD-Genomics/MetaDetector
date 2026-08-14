#!/usr/bin/env bash
set -Eeuo pipefail

MODE="docker"
PROFILE="test"
ROOT="/data/MD_Test_Dir"
IMAGE="metadetector:test"
THREADS="8"
MEMORY="16"
SRA_RUN="SRR22470068"
RUN_NAME="test_run"
YES="false"
MEGAN_MODE="auto"

MEGAN_MAP_URL="https://software-ab.cs.uni-tuebingen.de/download/megan7/megan-nr-r2.zip"
MEGAN_NUCL_URL="https://software-ab.cs.uni-tuebingen.de/download/megan7/megan-genome-r1.zip"

usage() {
cat <<'EOF'
MetaDetector fresh installer

Usage:
  ./install.sh [options]

Options:
  --mode docker                 Runtime mode (default: docker)
  --db-profile test|md          test = compact integration test resources
                                md   = full runtime/database layout only
  --root PATH                   Runtime root (default: /data/MD_Test_Dir)
  --image NAME                  Docker image (default: metadetector:test)
  --threads N                   Threads (default: 8)
  --memory GB                   Memory per stage in GB, integer only (default: 16)
  --sra-run SRR...              NCBI SRA run for test profile
                                (default: SRR22470068)
  --run-name NAME               Output run directory name
                                (default: test_run)
  --megan-mode auto|download|skip
                                auto: download for test, skip for md
                                download: fetch official MEGAN 7 mapping DBs
                                skip: do not provision MEGAN mapping DBs
  --yes                         Non-interactive
  -h, --help                    Show help

Notes:
- The Docker image must already exist.
- The image must contain corrected queue.sh/proc_assembly.sh, fasterq-dump,
  pigz, DIAMOND, MMseqs2, MEGAN CLI tools, and the local sbatch shim.
- The test profile provisions compact real NCBI references for BBMap,
  DIAMOND and MMseqs2.
- A full end-to-end test through MEGAN requires the official MEGAN protein
  and nucleotide mapping databases. These are substantially larger than
  the compact test references.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --db-profile) PROFILE="$2"; shift 2 ;;
        --root) ROOT="$2"; shift 2 ;;
        --image) IMAGE="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --memory) MEMORY="$2"; shift 2 ;;
        --sra-run) SRA_RUN="$2"; shift 2 ;;
        --run-name) RUN_NAME="$2"; shift 2 ;;
        --megan-mode) MEGAN_MODE="$2"; shift 2 ;;
        --yes) YES="true"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done

[[ "$MODE" == "docker" ]] || {
    echo "ERROR: only --mode docker is currently implemented" >&2
    exit 2
}

[[ "$PROFILE" == "test" || "$PROFILE" == "md" ]] || {
    echo "ERROR: --db-profile must be test or md" >&2
    exit 2
}

[[ "$THREADS" =~ ^[0-9]+$ ]] && (( THREADS > 0 )) || {
    echo "ERROR: --threads must be a positive integer" >&2
    exit 2
}

[[ "$MEMORY" =~ ^[0-9]+$ ]] && (( MEMORY > 0 )) || {
    echo "ERROR: --memory must be a positive integer in GB" >&2
    exit 2
}

[[ "$SRA_RUN" =~ ^[SED]RR[0-9]+$ ]] || {
    echo "ERROR: --sra-run must look like SRR123456, ERR..., or DRR..." >&2
    exit 2
}

[[ "$MEGAN_MODE" == "auto" || "$MEGAN_MODE" == "download" || "$MEGAN_MODE" == "skip" ]] || {
    echo "ERROR: --megan-mode must be auto, download, or skip" >&2
    exit 2
}

if [[ "$MEGAN_MODE" == "auto" ]]; then
    if [[ "$PROFILE" == "test" ]]; then
        MEGAN_MODE="download"
    else
        MEGAN_MODE="skip"
    fi
fi

command -v docker >/dev/null 2>&1 || {
    echo "ERROR: docker not found" >&2
    exit 1
}

docker image inspect "$IMAGE" >/dev/null 2>&1 || {
    echo "ERROR: Docker image not found: $IMAGE" >&2
    exit 1
}

echo "Checking Docker image..."
docker run --rm \
  --entrypoint /bin/bash \
  "$IMAGE" \
  -lc '
    set -e
    bash -n /opt/metadetector/queue.sh
    bash -n /opt/metadetector/proc_assembly.sh
    command -v fastqc >/dev/null
    command -v bbduk.sh >/dev/null
    command -v seqtk >/dev/null
    command -v pigz >/dev/null
    command -v diamond >/dev/null
    command -v mmseqs >/dev/null
    command -v sbatch >/dev/null
    command -v daa-meganizer >/dev/null || test -x /opt/megan/tools/daa-meganizer
    command -v blast2rma >/dev/null || test -x /opt/megan/tools/blast2rma
    test-sbatch-local >/dev/null
  '

if [[ "$PROFILE" == "test" ]]; then
    docker run --rm \
      --entrypoint /bin/bash \
      "$IMAGE" \
      -lc 'command -v fasterq-dump >/dev/null' || {
        echo "ERROR: image is missing fasterq-dump" >&2
        exit 1
    }
fi

mkdir -p \
  "$ROOT/config" \
  "$ROOT/input" \
  "$ROOT/output" \
  "$ROOT/databases" \
  "$ROOT/scratch/diamond" \
  "$ROOT/scratch/mmseqs" \
  "$ROOT/scripts" \
  "$ROOT/state/install-markers" \
  "$ROOT/state/install-logs" \
  "$ROOT/state/sbatch"

chmod -R u+rwX,g+rwX "$ROOT"

# Build container-visible passwd/group files for the invoking host identity.
CONTAINER_PASSWD="$ROOT/state/container_passwd"
CONTAINER_GROUP="$ROOT/state/container_group"

cp /etc/passwd "$CONTAINER_PASSWD"
cp /etc/group "$CONTAINER_GROUP"

HOST_UID="$(id -u)"
HOST_GID="$(id -g)"
HOST_USER="$(id -un 2>/dev/null || printf 'mduser')"

if ! awk -F: -v uid="$HOST_UID" '$3 == uid { found=1 } END { exit !found }' "$CONTAINER_PASSWD"; then
    if getent passwd "$HOST_UID" >/dev/null 2>&1; then
        getent passwd "$HOST_UID" >> "$CONTAINER_PASSWD"
    else
        printf '%s:x:%s:%s:MetaDetector Runtime User:/tmp:/bin/bash\n' \
          "$HOST_USER" "$HOST_UID" "$HOST_GID" >> "$CONTAINER_PASSWD"
    fi
fi

if ! awk -F: -v gid="$HOST_GID" '$3 == gid { found=1 } END { exit !found }' "$CONTAINER_GROUP"; then
    if getent group "$HOST_GID" >/dev/null 2>&1; then
        getent group "$HOST_GID" >> "$CONTAINER_GROUP"
    else
        printf 'mdhost:x:%s:\n' "$HOST_GID" >> "$CONTAINER_GROUP"
    fi
fi

chmod 0644 "$CONTAINER_PASSWD" "$CONTAINER_GROUP"

write_md_config() {
cat > "$ROOT/config/md.config" <<'EOF'
# Generated MetaDetector Docker runtime config

md_partition="local"

md_dbdir="/database"
md_meganpath="/opt/megan/tools"

md_diamond_args="--block-size 2 --index-chunks 1 --tmpdir /scratch/diamond"
md_diamond_dbdir="/database/diamond/nr.dmnd"

md_mmseqs_dbdir="/database/mmseqs/core_nt/core_nt"
md_mmseqs_tmpdir="/scratch/mmseqs"

# queue.sh uses md_datadir to locate proc_assembly.sh
md_datadir="/opt/metadetector"

md_docker_execdir="/opt/metadetector"
md_docker_tmpdir="/scratch"
md_docker_megandir="/opt/megan"
md_docker_slurmconf=""

md_human_ref="/database/bbmap/human/human.fa"
md_silva_ref="/database/bbmap/silva/silva.fa"
md_contam_ref="/database/bbmap/contam/contam.fa"
EOF
}

write_run_env() {
cat > "$ROOT/config/run.env" <<EOF
MD_ROOT="$ROOT"
MD_IMAGE="$IMAGE"
MD_PROFILE="$PROFILE"
MD_THREADS="$THREADS"
MD_MEMORY="$MEMORY"
MD_RUN_NAME="$RUN_NAME"
MD_ASSEMBLY_FLAG="m"
MD_ASSEMBLY_TYPE="s"
MD_SRA_RUN="$SRA_RUN"
MD_MEGAN_MODE="$MEGAN_MODE"
EOF

    if [[ "$PROFILE" == "test" ]]; then
        echo 'MD_STAGES="1,2,3a,4,5"' >> "$ROOT/config/run.env"
        echo 'MD_EXPECTED_STAGES="1 2 3a 4 5 6 9 12 14 16 17 19 21 22 25"' >> "$ROOT/config/run.env"
    else
        echo 'MD_STAGES=""' >> "$ROOT/config/run.env"
        echo 'MD_EXPECTED_STAGES=""' >> "$ROOT/config/run.env"
    fi
}

download_test_reads() {
    local marker="$ROOT/state/install-markers/test_reads_${SRA_RUN}.done"
    local log="$ROOT/state/install-logs/test_reads_${SRA_RUN}.log"

    if [[ -e "$marker" \
       && -s "$ROOT/input/${SRA_RUN}_R1.fastq.gz" \
       && -s "$ROOT/input/${SRA_RUN}_R2.fastq.gz" ]]; then
        echo "SKIP NCBI test reads: $SRA_RUN"
        return
    fi

    echo "Downloading NCBI test reads: $SRA_RUN"

    rm -f \
      "$ROOT/input/${SRA_RUN}_R1.fastq.gz" \
      "$ROOT/input/${SRA_RUN}_R2.fastq.gz"

    rm -rf "$ROOT/input/.sra_tmp_${SRA_RUN}"
    mkdir -p "$ROOT/input/.sra_tmp_${SRA_RUN}"

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      -v "$ROOT/state/container_passwd:/etc/passwd:ro" \
      -v "$ROOT/state/container_group:/etc/group:ro" \
      -e HOME=/tmp \
      -e TMPDIR=/scratch \
      -v "$ROOT/input:/input" \
      -v "$ROOT/scratch:/scratch" \
      --entrypoint /bin/bash \
      "$IMAGE" \
      -lc '
        set -Eeuo pipefail

        run="'"$SRA_RUN"'"
        tmp="/input/.sra_tmp_${run}"

        fasterq-dump \
          "$run" \
          --split-files \
          --threads "'"$THREADS"'" \
          --outdir "$tmp" \
          --temp /scratch

        test -s "$tmp/${run}_1.fastq"
        test -s "$tmp/${run}_2.fastq"

        awk "NR % 4 == 3 {\$0=\"+\"} {print}" \
          "$tmp/${run}_1.fastq" \
          | pigz -p "'"$THREADS"'" \
          > "/input/${run}_R1.fastq.gz"

        awk "NR % 4 == 3 {\$0=\"+\"} {print}" \
          "$tmp/${run}_2.fastq" \
          | pigz -p "'"$THREADS"'" \
          > "/input/${run}_R2.fastq.gz"

        test -s "/input/${run}_R1.fastq.gz"
        test -s "/input/${run}_R2.fastq.gz"

        zcat "/input/${run}_R1.fastq.gz" | sed -n "3p" | grep -qx "+"
        zcat "/input/${run}_R2.fastq.gz" | sed -n "3p" | grep -qx "+"

        rm -rf "$tmp"
      ' 2>&1 | tee "$log"

    test -s "$ROOT/input/${SRA_RUN}_R1.fastq.gz"
    test -s "$ROOT/input/${SRA_RUN}_R2.fastq.gz"
    gzip -t "$ROOT/input/${SRA_RUN}_R1.fastq.gz"
    gzip -t "$ROOT/input/${SRA_RUN}_R2.fastq.gz"

    printf '/input/%s_R1.fastq.gz\n/input/%s_R2.fastq.gz\n' \
      "$SRA_RUN" "$SRA_RUN" \
      > "$ROOT/config/reads.txt"

    printf '%s\n' "$SRA_RUN" > "$ROOT/config/sample_name.txt"

    touch "$marker"
}

provision_test_databases() {
    local marker="$ROOT/state/install-markers/test_databases.done"
    local log="$ROOT/state/install-logs/test_databases.log"

    local human="$ROOT/databases/bbmap/human/human.fa"
    local contam="$ROOT/databases/bbmap/contam/contam.fa"
    local silva="$ROOT/databases/bbmap/silva/silva.fa"

    local diamond_faa="$ROOT/databases/diamond/nr_test.faa"
    local diamond_db="$ROOT/databases/diamond/nr.dmnd"
    local positive_query="$ROOT/databases/diamond/positive_query.fa"
    local diamond_smoke="$ROOT/databases/diamond/diamond_smoke.tsv"

    local mmseqs_fa="$ROOT/databases/mmseqs/core_nt/core_nt_test.fa"
    local mmseqs_db="$ROOT/databases/mmseqs/core_nt/core_nt"
    local mmseqs_smoke="$ROOT/databases/mmseqs/core_nt/mmseqs_smoke.tsv"

    local megan_map="$ROOT/databases/megan/megan-map.db"
    local megan_nucl="$ROOT/databases/megan/megan-nucl.db"

    local complete=true
    for f in \
      "$human" "$contam" "$silva" \
      "$diamond_faa" "$diamond_db" "$positive_query" "$diamond_smoke" \
      "$mmseqs_fa" "$mmseqs_db" "$mmseqs_db.dbtype" "$mmseqs_smoke"
    do
        [[ -s "$f" ]] || complete=false
    done

    if [[ "$MEGAN_MODE" == "download" ]]; then
        [[ -s "$megan_map" ]] || complete=false
        [[ -s "$megan_nucl" ]] || complete=false
    fi

    if [[ -e "$marker" && "$complete" == "true" ]]; then
        echo "SKIP compact test databases"
        return
    fi

    echo "Provisioning compact sample-derived test databases..."

    rm -f "$marker"
    mkdir -p \
      "$ROOT/databases/bbmap/human" \
      "$ROOT/databases/bbmap/contam" \
      "$ROOT/databases/bbmap/silva" \
      "$ROOT/databases/diamond" \
      "$ROOT/databases/mmseqs/core_nt" \
      "$ROOT/databases/megan"

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      -v "$ROOT/state/container_passwd:/etc/passwd:ro" \
      -v "$ROOT/state/container_group:/etc/group:ro" \
      -e HOME=/tmp \
      -e TEST_RUN="$SRA_RUN" \
      -v "$ROOT/input:/input:ro" \
      -v "$ROOT/databases:/database" \
      --entrypoint /bin/bash \
      "$IMAGE" \
      -lc '
        set -Eeuo pipefail

        python - <<'"'"'PY'"'"'
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
import gzip
import os
import time

MAX_RETRIES = 6
run = os.environ["TEST_RUN"]

def download(url, dst, min_size=100):
    dst = Path(dst)

    if dst.exists() and dst.stat().st_size >= min_size:
        print(f"SKIP existing {dst}\t{dst.stat().st_size} bytes")
        return dst.read_bytes()

    dst.parent.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with urlopen(url, timeout=120) as r:
                data = r.read()

            if len(data) < min_size:
                raise RuntimeError(
                    f"Downloaded file too small for {dst}: {len(data)} bytes"
                )

            dst.write_bytes(data)
            print(f"{dst}\t{dst.stat().st_size} bytes")
            time.sleep(1.0)
            return data

        except HTTPError as e:
            if e.code == 429 and attempt < MAX_RETRIES:
                retry_after = e.headers.get("Retry-After")
                try:
                    delay = int(retry_after) if retry_after else min(60, 5 * attempt)
                except ValueError:
                    delay = min(60, 5 * attempt)

                print(
                    f"NCBI rate limit for {dst}; "
                    f"retry {attempt}/{MAX_RETRIES} in {delay}s"
                )
                time.sleep(delay)
                continue
            raise

        except URLError:
            if attempt < MAX_RETRIES:
                delay = min(60, 5 * attempt)
                print(
                    f"NCBI connection error for {dst}; "
                    f"retry {attempt}/{MAX_RETRIES} in {delay}s"
                )
                time.sleep(delay)
                continue
            raise

refs = {
    "/database/bbmap/human/human.fa":
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        "?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text",

    "/database/bbmap/contam/contam.fa":
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        "?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text",

    "/database/bbmap/silva/silva.fa":
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        "?db=nuccore&id=NR_024570.1&rettype=fasta&retmode=text",
}

for dst, url in refs.items():
    data = download(url, dst)
    if not data.startswith(b">"):
        raise RuntimeError(f"NCBI did not return FASTA for {dst}")

r1 = Path(f"/input/{run}_R1.fastq.gz")
r2 = Path(f"/input/{run}_R2.fastq.gz")

for p in (r1, r2):
    if not p.exists() or p.stat().st_size == 0:
        raise RuntimeError(f"Missing test FASTQ: {p}")

mmseqs_fa = Path("/database/mmseqs/core_nt/core_nt_test.fa")
diamond_faa = Path("/database/diamond/nr_test.faa")
positive_query = Path("/database/diamond/positive_query.fa")
manifest = Path("/database/diamond/TEST_ONLY_README.txt")

mmseqs_fa.parent.mkdir(parents=True, exist_ok=True)
diamond_faa.parent.mkdir(parents=True, exist_ok=True)

comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")

codon_table = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"X","TAG":"X",
    "TGT":"C","TGC":"C","TGA":"X","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def translate(seq, frame):
    seq = seq.upper()
    aa = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa.append(codon_table.get(codon, "X"))
    return "".join(aa)

def revcomp(seq):
    return seq.translate(comp)[::-1]

def fastq_records(path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                raise RuntimeError(f"Truncated FASTQ: {path}")
            yield h[1:].strip().split()[0], seq

n_reads = 0
n_proteins = 0
positive_written = 0

with mmseqs_fa.open("w") as nuc, diamond_faa.open("w") as prot, positive_query.open("w") as pq:
    for mate_no, path in enumerate((r1, r2), start=1):
        for read_id, seq in fastq_records(path):
            n_reads += 1
            clean_id = f"mdtest_{mate_no}_{n_reads}"

            nuc.write(f">{clean_id}\n{seq}\n")

            if positive_written < 200:
                pq.write(f">{clean_id}\n{seq}\n")
                positive_written += 1

            strands = (("f", seq), ("r", revcomp(seq)))
            for strand_name, strand_seq in strands:
                for frame in range(3):
                    aa = translate(strand_seq, frame)
                    if len(aa) < 20:
                        continue

                    n_proteins += 1

                    # ref|...| allows MEGAN accession-tag parsing while the
                    # trailing fields keep DIAMOND subject IDs unique.
                    prot.write(
                        f">ref|WP_013615770.1|mdtest_{n_proteins}_{strand_name}{frame+1}\n"
                        f"{aa}\n"
                    )

if n_reads == 0 or n_proteins == 0 or positive_written == 0:
    raise RuntimeError(
        f"Failed building sample-derived test DBs: "
        f"reads={n_reads}, proteins={n_proteins}, positive={positive_written}"
    )

manifest.write_text(
    "MetaDetector TEST-ONLY DIAMOND database\n"
    "=======================================\n"
    "Purpose: installation/integration smoke testing only.\n"
    "Protein sequences are six-frame translations of the test FASTQs.\n"
    "The RefSeq-style accession token in subject headers is used only to\n"
    "exercise MEGAN accession-mapping code paths. These sequences MUST NOT\n"
    "be interpreted biologically as WP_013615770.1.\n"
    f"Source run: {run}\n"
    f"Nucleotide reads in MMseqs FASTA: {n_reads}\n"
    f"Protein entries in DIAMOND FASTA: {n_proteins}\n"
)

print(f"{mmseqs_fa}\t{mmseqs_fa.stat().st_size} bytes\t{n_reads} sequences")
print(f"{diamond_faa}\t{diamond_faa.stat().st_size} bytes\t{n_proteins} proteins")
print(f"{positive_query}\t{positive_query.stat().st_size} bytes\t{positive_written} sequences")
PY

        for f in \
          /database/bbmap/human/human.fa \
          /database/bbmap/contam/contam.fa \
          /database/bbmap/silva/silva.fa \
          /database/diamond/nr_test.faa \
          /database/diamond/positive_query.fa \
          /database/mmseqs/core_nt/core_nt_test.fa
        do
            test -s "$f"
            head -1 "$f" | grep -q "^>"
        done

        rm -f /database/diamond/nr.dmnd
        diamond makedb \
          --in /database/diamond/nr_test.faa \
          --db /database/diamond/nr
        test -s /database/diamond/nr.dmnd

        rm -f \
          /database/mmseqs/core_nt/core_nt \
          /database/mmseqs/core_nt/core_nt.dbtype \
          /database/mmseqs/core_nt/core_nt.index \
          /database/mmseqs/core_nt/core_nt.lookup \
          /database/mmseqs/core_nt/core_nt.source

        mmseqs createdb \
          /database/mmseqs/core_nt/core_nt_test.fa \
          /database/mmseqs/core_nt/core_nt

        test -s /database/mmseqs/core_nt/core_nt
        test -s /database/mmseqs/core_nt/core_nt.dbtype

        rm -f /database/diamond/diamond_smoke.tsv
        diamond blastx \
          --query /database/diamond/positive_query.fa \
          --db /database/diamond/nr.dmnd \
          --out /database/diamond/diamond_smoke.tsv \
          --outfmt 6 \
          --max-target-seqs 1 \
          --threads 4

        test -s /database/diamond/diamond_smoke.tsv

        rm -f /database/mmseqs/core_nt/mmseqs_smoke.tsv
        rm -rf /tmp/md_mmseqs_smoke

        mmseqs easy-search \
          /database/diamond/positive_query.fa \
          /database/mmseqs/core_nt/core_nt \
          /database/mmseqs/core_nt/mmseqs_smoke.tsv \
          /tmp/md_mmseqs_smoke \
          --search-type 3 \
          -e 1e-20 \
          --threads 4

        test -s /database/mmseqs/core_nt/mmseqs_smoke.tsv

        echo "DIAMOND positive-control hits: $(wc -l < /database/diamond/diamond_smoke.tsv)"
        echo "MMseqs positive-control hits:  $(wc -l < /database/mmseqs/core_nt/mmseqs_smoke.tsv)"
      ' 2>&1 | tee "$log"

    for f in \
      "$human" "$contam" "$silva" \
      "$diamond_faa" "$diamond_db" "$positive_query" "$diamond_smoke" \
      "$mmseqs_fa" "$mmseqs_db" "$mmseqs_db.dbtype" "$mmseqs_smoke"
    do
        test -s "$f"
    done

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      -v "$ROOT/state/container_passwd:/etc/passwd:ro" \
      -v "$ROOT/state/container_group:/etc/group:ro" \
      -v "$ROOT/databases:/database:ro" \
      --entrypoint /bin/bash \
      "$IMAGE" \
      -lc '
        set -e
        diamond dbinfo --db /database/diamond/nr.dmnd >/dev/null
        mmseqs dbtype /database/mmseqs/core_nt/core_nt >/dev/null
        test -s /database/diamond/diamond_smoke.tsv
        test -s /database/mmseqs/core_nt/mmseqs_smoke.tsv
      '

    if [[ "$MEGAN_MODE" == "download" ]]; then
        provision_megan_databases
    else
        echo "MEGAN provisioning skipped."
    fi

    touch "$marker"
}


provision_megan_databases() {
    local marker="$ROOT/state/install-markers/megan_databases.done"
    local log="$ROOT/state/install-logs/megan_databases.log"

    local megan_map="$ROOT/databases/megan/megan-map.db"
    local megan_nucl="$ROOT/databases/megan/megan-nucl.db"

    if [[ -e "$marker" && -s "$megan_map" && -s "$megan_nucl" ]]; then
        echo "SKIP MEGAN mapping databases"
        return
    fi

    echo
    echo "Provisioning official MEGAN mapping databases."
    echo "NOTE: these downloads are much larger than the compact test references."

    rm -f "$marker"
    mkdir -p "$ROOT/databases/megan"

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      -v "$ROOT/state/container_passwd:/etc/passwd:ro" \
      -v "$ROOT/state/container_group:/etc/group:ro" \
      -e HOME=/tmp \
      -v "$ROOT/databases:/database" \
      --entrypoint /bin/bash \
      "$IMAGE" \
      -lc '
        set -Eeuo pipefail

        MAP_URL="'"$MEGAN_MAP_URL"'"
        NUCL_URL="'"$MEGAN_NUCL_URL"'"

        python - <<'"'"'PY'"'"'
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile
import shutil

jobs = [
    (
        "'"$MEGAN_MAP_URL"'",
        Path("/database/megan/megan-nr-r2.zip"),
        "megan-nr-r2.mdb",
        Path("/database/megan/megan-map.db"),
    ),
    (
        "'"$MEGAN_NUCL_URL"'",
        Path("/database/megan/megan-genome-r1.zip"),
        "megan-genome-r1.mdb",
        Path("/database/megan/megan-nucl.db"),
    ),
]

for url, zip_path, expected_member, final_path in jobs:
    if final_path.exists() and final_path.stat().st_size > 0:
        print(f"SKIP existing {final_path}")
        continue

    print(f"Downloading {url}")
    with urlopen(url, timeout=300) as response, zip_path.open("wb") as out:
        shutil.copyfileobj(response, out, length=16 * 1024 * 1024)

    if zip_path.stat().st_size == 0:
        raise RuntimeError(f"Empty download: {zip_path}")

    with ZipFile(zip_path) as z:
        members = z.namelist()

        member = None
        for name in members:
            if name == expected_member or name.endswith("/" + expected_member):
                member = name
                break

        if member is None:
            mdb_members = [n for n in members if n.endswith(".mdb")]
            if len(mdb_members) != 1:
                raise RuntimeError(
                    f"Expected {expected_member} in {zip_path}; found: {mdb_members}"
                )
            member = mdb_members[0]

        with z.open(member) as src, final_path.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=16 * 1024 * 1024)

    if final_path.stat().st_size == 0:
        raise RuntimeError(f"Empty extracted database: {final_path}")

    zip_path.unlink(missing_ok=True)
    print(f"{final_path}\t{final_path.stat().st_size} bytes")
PY

        test -s /database/megan/megan-map.db
        test -s /database/megan/megan-nucl.db
      ' 2>&1 | tee "$log"

    test -s "$megan_map"
    test -s "$megan_nucl"

    touch "$marker"
}

write_runner() {
cat > "$ROOT/scripts/run_metadetector.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

ENV_FILE="${1:-}"
[[ -n "$ENV_FILE" && -f "$ENV_FILE" ]] || {
    echo "Usage: $0 /path/to/run.env" >&2
    exit 2
}

# shellcheck disable=SC1090
source "$ENV_FILE"

: "${MD_ROOT:?}"
: "${MD_IMAGE:?}"
: "${MD_THREADS:?}"
: "${MD_MEMORY:?}"
: "${MD_RUN_NAME:?}"
: "${MD_ASSEMBLY_FLAG:?}"
: "${MD_ASSEMBLY_TYPE:?}"

RUN_DIR="$MD_ROOT/output/$MD_RUN_NAME"
LOG_DIR="$RUN_DIR/logs"
MASTER_LOG="$LOG_DIR/run.log"
META_FILE="$RUN_DIR/run_metadata.txt"
PSEUDO_SLURM_DIR="$MD_ROOT/state/sbatch"
SAMPLE="${MD_SRA_RUN:-}"
[[ -z "$SAMPLE" && -f "$MD_ROOT/config/sample_name.txt" ]] && SAMPLE="$(cat "$MD_ROOT/config/sample_name.txt")"
SAMPLE_DIR="$RUN_DIR/$SAMPLE"
STATUS_DIR="$SAMPLE_DIR/status_log"

mkdir -p "$RUN_DIR" "$LOG_DIR" "$PSEUDO_SLURM_DIR" \
  "$MD_ROOT/scratch/diamond" "$MD_ROOT/scratch/mmseqs"

rm -f "$RUN_DIR/run.finished" "$RUN_DIR/run.failed" "$RUN_DIR/run.interrupted"
rm -f "$MASTER_LOG"

stage_args=()
[[ -n "${MD_STAGES:-}" ]] && stage_args=(-s "$MD_STAGES")

stage_label() {
    case "$1" in
        1)   echo "Input preparation" ;;
        2)   echo "Read QC / trimming" ;;
        3a)  echo "Host removal" ;;
        4)   echo "Contaminant removal" ;;
        5)   echo "rRNA removal" ;;
        6)   echo "Assembly" ;;
        9)   echo "Assembly processing" ;;
        12)  echo "DIAMOND contigs" ;;
        14)  echo "DIAMOND reads" ;;
        16)  echo "MMseqs nucleotide search" ;;
        17)  echo "MEGAN contigs" ;;
        19)  echo "MEGAN reads" ;;
        21)  echo "MEGAN nucleotide" ;;
        22)  echo "MetaQUAST" ;;
        25)  echo "Final QC" ;;
        *)   echo "Pipeline stage" ;;
    esac
}

stamp() { date '+%H:%M:%S'; }

{
    echo "MetaDetector run"
    echo "started=$(date -Is)"
    echo "run_name=$MD_RUN_NAME"
    echo "sample=$SAMPLE"
    echo "profile=${MD_PROFILE:-unknown}"
    echo "image=$MD_IMAGE"
    echo "threads=$MD_THREADS"
    echo "memory=$MD_MEMORY"
    echo "stages=${MD_STAGES:-all}"
    echo "expected_stages=${MD_EXPECTED_STAGES:-unknown}"
} > "$META_FILE"

cat <<HDR
============================================================
 MetaDetector Test
============================================================
Run:       $MD_RUN_NAME
Sample:    ${SAMPLE:-unknown}
Image:     $MD_IMAGE
Threads:   $MD_THREADS
Memory:    $MD_MEMORY GB
Results:   $RUN_DIR
Log:       $MASTER_LOG
============================================================
Starting MetaDetector...
HDR

CONTAINER_NAME="metadetector_${MD_RUN_NAME}_$$"
cleanup_container() { docker rm -f "$CONTAINER_NAME" >/dev/null 2>&1 || true; }
interrupt_run() {
    echo
    echo "[$(stamp)] INTERRUPTED - stopping MetaDetector container"
    cleanup_container
    touch "$RUN_DIR/run.interrupted"
    exit 130
}
trap interrupt_run INT TERM TSTP

rm -rf "$PSEUDO_SLURM_DIR"
mkdir -p "$PSEUDO_SLURM_DIR"

set +e
docker run \
  --name "$CONTAINER_NAME" \
  --user "$(id -u):$(id -g)" \
  -v "$MD_ROOT/state/container_passwd:/etc/passwd:ro" \
  -v "$MD_ROOT/state/container_group:/etc/group:ro" \
  -v "$MD_ROOT/config/md.config:/opt/metadetector/md.config:ro" \
  -v "$MD_ROOT/config/reads.txt:/work/reads.txt:ro" \
  -v "$MD_ROOT/input:/input:ro" \
  -v "$MD_ROOT/output:/output" \
  -v "$MD_ROOT/databases:/database:ro" \
  -v "$MD_ROOT/scratch:/scratch" \
  -v "$MD_ROOT/state:/state" \
  -e HOME=/tmp \
  -e MD_LOCAL_SBATCH_STATE_DIR="/state/sbatch" \
  -e TMPDIR=/scratch \
  --entrypoint /bin/bash \
  "$MD_IMAGE" \
  -lc '
    set -Eeuo pipefail
    cd /opt/metadetector
    ./queue.sh \
      -r /work/reads.txt \
      -o "/output/'"$MD_RUN_NAME"'" \
      -a "'"$MD_ASSEMBLY_FLAG"'" \
      -j "'"$MD_ASSEMBLY_TYPE"'" \
      -m "'"$MD_MEMORY"'" \
      -t "'"$MD_THREADS"'" \
      '"$(printf '%q ' "${stage_args[@]}")"'
  ' > >(tee -a "$MASTER_LOG") 2>&1 &
docker_pid=$!
set -e

declare -A shown=()
while kill -0 "$docker_pid" 2>/dev/null; do
    for stage in ${MD_EXPECTED_STAGES:-}; do
        if [[ -z "${shown[$stage]:-}" && -e "$STATUS_DIR/stage${stage}.finished" ]]; then
            printf '[%s] PASS  Stage %-3s %s\n' "$(stamp)" "$stage" "$(stage_label "$stage")"
            shown[$stage]=1
        elif [[ -z "${shown[$stage]:-}" && -e "$STATUS_DIR/stage${stage}.not_finished" ]]; then
            printf '[%s] FAIL  Stage %-3s %s\n' "$(stamp)" "$stage" "$(stage_label "$stage")"
            shown[$stage]=fail
        fi
    done
    sleep 2
done

set +e
wait "$docker_pid"
rc=$?
set -e
cleanup_container
trap - INT TERM TSTP

for stage in ${MD_EXPECTED_STAGES:-}; do
    if [[ -z "${shown[$stage]:-}" && -e "$STATUS_DIR/stage${stage}.finished" ]]; then
        printf '[%s] PASS  Stage %-3s %s\n' "$(stamp)" "$stage" "$(stage_label "$stage")"
        shown[$stage]=1
    fi
done

missing=()
for stage in ${MD_EXPECTED_STAGES:-}; do
    [[ -e "$STATUS_DIR/stage${stage}.finished" ]] || missing+=("$stage")
done

completion_ok=0
if [[ -s "$STATUS_DIR/assembly_pipeline.finished" ]] && \
   grep -q 'assembly pipeline completed' "$STATUS_DIR/assembly_pipeline.finished"; then
    completion_ok=1
fi

scheduler_ok=1
scheduler_bad=()
shopt -s nullglob
rc_files=("$PSEUDO_SLURM_DIR"/*.rc)
shopt -u nullglob
if (( ${#rc_files[@]} == 0 )); then
    scheduler_ok=0
else
    for f in "${rc_files[@]}"; do
        value="$(tr -d '[:space:]' < "$f")"
        if [[ "$value" != "0" ]]; then
            scheduler_ok=0
            scheduler_bad+=("$(basename "$f" .rc):$value")
        fi
    done
fi

daa_ok=1
daa_summary="not checked"
if [[ "${MD_PROFILE:-}" == "test" && -n "$SAMPLE" ]]; then
    contigs_daa="$SAMPLE_DIR/blast/${SAMPLE}_metaspades_contigs_blastx.daa"
    reads_daa="$SAMPLE_DIR/blast/${SAMPLE}_reads_blastx.daa"
    if [[ -s "$contigs_daa" && -s "$reads_daa" ]]; then
        set +e
        daa_summary="$(docker run --rm \
          --user "$(id -u):$(id -g)" \
          -v "$MD_ROOT/output:/output:ro" \
          --entrypoint /bin/bash "$MD_IMAGE" -lc '
            set -e
            c="/output/'"$MD_RUN_NAME/$SAMPLE"'/blast/'"$SAMPLE"'_metaspades_contigs_blastx.daa"
            r="/output/'"$MD_RUN_NAME/$SAMPLE"'/blast/'"$SAMPLE"'_reads_blastx.daa"
            diamond view --daa "$c" --out /tmp/c.tsv >/dev/null 2>&1
            diamond view --daa "$r" --out /tmp/r.tsv >/dev/null 2>&1
            c_n=$(wc -l < /tmp/c.tsv)
            r_n=$(wc -l < /tmp/r.tsv)
            test "$c_n" -gt 0
            test "$r_n" -gt 0
            printf "contigs=%s reads=%s" "$c_n" "$r_n"
          ' 2>/dev/null)"
        daa_rc=$?
        set -e
        [[ $daa_rc -eq 0 ]] || daa_ok=0
    else
        daa_ok=0
    fi
fi

if [[ $rc -ne 0 || ${#missing[@]} -gt 0 || $completion_ok -ne 1 || $scheduler_ok -ne 1 || $daa_ok -ne 1 ]]; then
    rc=1
fi

{
    echo "finished=$(date -Is)"
    echo "exit_code=$rc"
    echo "completion_marker=$completion_ok"
    echo "scheduler_ok=$scheduler_ok"
    echo "daa_ok=$daa_ok"
    echo "daa_summary=$daa_summary"
} >> "$META_FILE"

echo
echo "============================================================"
if [[ $rc -eq 0 ]]; then
    touch "$RUN_DIR/run.finished"
    rm -f "$RUN_DIR/run.failed"
    echo " MetaDetector test: PASS"
else
    touch "$RUN_DIR/run.failed"
    rm -f "$RUN_DIR/run.finished"
    echo " MetaDetector test: FAIL"
fi
echo "============================================================"
printf 'Required stages:          %s/%s PASS\n' \
  "$(( $(wc -w <<<"${MD_EXPECTED_STAGES:-}") - ${#missing[@]} ))" \
  "$(wc -w <<<"${MD_EXPECTED_STAGES:-}")"
[[ $completion_ok -eq 1 ]] && echo "Final completion marker:   PASS" || echo "Final completion marker:   FAIL"
[[ $scheduler_ok -eq 1 ]] && echo "Local scheduler jobs:      PASS" || echo "Local scheduler jobs:      FAIL"
if [[ "${MD_PROFILE:-}" == "test" ]]; then
    [[ $daa_ok -eq 1 ]] && echo "DIAMOND DAA validation:    PASS ($daa_summary)" || echo "DIAMOND DAA validation:    FAIL"
fi

if (( ${#missing[@]} > 0 )); then
    echo "Missing stages:            ${missing[*]}"
fi
if (( ${#scheduler_bad[@]} > 0 )); then
    echo "Failed scheduler jobs:     ${scheduler_bad[*]}"
fi

echo "Results:                   $RUN_DIR"
echo "Detailed log:              $MASTER_LOG"

echo
if [[ $rc -eq 0 ]]; then
    echo "Test completed successfully. No additional debugging commands are required."
else
    echo "Failure details:"
    grep -HniE 'fatal|traceback|exception|command not found|no such file|out of memory|oom|sequence set index out of bounds' \
      "$SAMPLE_DIR"/log/stage*.err "$SAMPLE_DIR"/log/stage*.out 2>/dev/null | tail -30 || true
fi

exit "$rc"
EOF

chmod 0755 "$ROOT/scripts/run_metadetector.sh"

cat > "$ROOT/run_test.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$ROOT/scripts/run_metadetector.sh" "$ROOT/config/run.env"
EOF
chmod 0755 "$ROOT/run_test.sh"
}

write_validator() {
cat > "$ROOT/scripts/validate_install.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# shellcheck disable=SC1091
source "$ROOT/config/run.env"

fail=0

ok() {
    printf 'PASS: %s\n' "$*"
}

bad() {
    printf 'FAIL: %s\n' "$*" >&2
    fail=1
}

check_file() {
    local f="$1"
    [[ -s "$f" ]] && ok "$f" || bad "$f"
}

echo "===== MetaDetector installation validation ====="

[[ -f "$ROOT/config/md.config" ]] && ok "md.config" || bad "md.config"
[[ -f "$ROOT/config/run.env" ]] && ok "run.env" || bad "run.env"
[[ -x "$ROOT/scripts/run_metadetector.sh" ]] && ok "run script" || bad "run script"

if docker image inspect "$MD_IMAGE" >/dev/null 2>&1; then
    ok "Docker image: $MD_IMAGE"
else
    bad "Docker image: $MD_IMAGE"
fi

if docker run --rm \
    --entrypoint /bin/bash \
    "$MD_IMAGE" \
    -lc '
      set -e
      bash -n /opt/metadetector/queue.sh
      bash -n /opt/metadetector/proc_assembly.sh
      command -v fastqc >/dev/null
      command -v bbduk.sh >/dev/null
      command -v seqtk >/dev/null
      command -v pigz >/dev/null
      command -v diamond >/dev/null
      command -v mmseqs >/dev/null
      command -v sbatch >/dev/null
      test-sbatch-local
    '
then
    ok "container smoke test"
else
    bad "container smoke test"
fi

for f in \
  "$ROOT/databases/bbmap/human/human.fa" \
  "$ROOT/databases/bbmap/contam/contam.fa" \
  "$ROOT/databases/bbmap/silva/silva.fa" \
  "$ROOT/databases/diamond/nr.dmnd" \
  "$ROOT/databases/diamond/diamond_smoke.tsv" \
  "$ROOT/databases/mmseqs/core_nt/core_nt" \
  "$ROOT/databases/mmseqs/core_nt/core_nt.dbtype" \
  "$ROOT/databases/mmseqs/core_nt/mmseqs_smoke.tsv"
do
    check_file "$f"
done

if [[ "${MD_MEGAN_MODE:-skip}" == "download" || "$MD_PROFILE" == "md" ]]; then
    check_file "$ROOT/databases/megan/megan-map.db"
    check_file "$ROOT/databases/megan/megan-nucl.db"
fi

if [[ "$MD_PROFILE" == "test" ]]; then
    check_file "$ROOT/input/${MD_SRA_RUN}_R1.fastq.gz"
    check_file "$ROOT/input/${MD_SRA_RUN}_R2.fastq.gz"
fi

exit "$fail"
EOF

chmod 0755 "$ROOT/scripts/validate_install.sh"
}

write_stage_reporter() {
cat > "$ROOT/scripts/report_stages.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/config/run.env"
sample="${MD_SRA_RUN:-$(cat "$ROOT/config/sample_name.txt" 2>/dev/null || true)}"
run="$ROOT/output/$MD_RUN_NAME/$sample"
status="$run/status_log"

pass=0
total=0
echo "===== METADETECTOR STATUS ====="
for s in ${MD_EXPECTED_STAGES:-}; do
    total=$((total+1))
    if [[ -e "$status/stage${s}.finished" ]]; then
        printf "stage%-5s PASS\n" "$s"
        pass=$((pass+1))
    elif [[ -e "$status/stage${s}.not_finished" ]]; then
        printf "stage%-5s FAIL\n" "$s"
    else
        printf "stage%-5s WAITING\n" "$s"
    fi
done

echo
printf "Stages: %s/%s complete\n" "$pass" "$total"
if [[ -s "$status/assembly_pipeline.finished" ]] && grep -q 'assembly pipeline completed' "$status/assembly_pipeline.finished"; then
    echo "Final completion: PASS"
else
    echo "Final completion: WAITING"
fi

if [[ -e "$ROOT/output/$MD_RUN_NAME/run.finished" ]]; then
    echo "Test: PASS"
elif [[ -e "$ROOT/output/$MD_RUN_NAME/run.failed" ]]; then
    echo "Test: FAIL"
else
    echo "Test: RUNNING"
fi
EOF
chmod 0755 "$ROOT/scripts/report_stages.sh"
}

write_md_config
write_run_env
write_runner
write_validator
write_stage_reporter

if [[ "$PROFILE" == "test" ]]; then
    download_test_reads
    provision_test_databases
else
    : > "$ROOT/config/reads.txt"
    mkdir -p \
      "$ROOT/databases/bbmap/human" \
      "$ROOT/databases/bbmap/contam" \
      "$ROOT/databases/bbmap/silva" \
      "$ROOT/databases/diamond" \
      "$ROOT/databases/mmseqs/core_nt" \
      "$ROOT/databases/megan"
fi

echo
echo "============================================================"
echo "MetaDetector installation complete"
echo "============================================================"
echo "Root       : $ROOT"
echo "Profile    : $PROFILE"
echo "Image      : $IMAGE"
echo "Run name   : $RUN_NAME"
echo "MEGAN mode : $MEGAN_MODE"

if [[ "$PROFILE" == "test" ]]; then
    echo "NCBI run   : $SRA_RUN"
fi

echo
echo "Validate:"
echo "  $ROOT/scripts/validate_install.sh"

if [[ "$PROFILE" == "test" ]]; then
    echo
    echo "Run:"
    echo "  $ROOT/run_test.sh"
    echo
    echo "The runner prints live stage progress and performs final validation automatically."
fi
