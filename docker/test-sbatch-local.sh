#!/usr/bin/env bash
set -Eeuo pipefail

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

export MD_LOCAL_SBATCH_STATE_DIR="$tmp/state"

cat > "$tmp/job.sh" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=md-local-test
#SBATCH --output=$tmp/job.out
#SBATCH --error=$tmp/job.err
#SBATCH --open-mode=append
echo LOCAL_SBATCH_OK
EOF

out="$(sbatch "$tmp/job.sh")"
grep -q '^Submitted batch job ' <<<"$out"
grep -q 'LOCAL_SBATCH_OK' "$tmp/job.out"
echo "local sbatch self-test: PASS"
