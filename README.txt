MetaDetector installer v4

Changes:
- removes queue.sh R1/R2 symlink creation for short-read runs
- resolves mates only after the current file is copied
- mounts a runtime-patched queue.sh into the container
- stores pseudo-SLURM state under output/<run>/pseudo_slurm/

Run:
  rm -rf /data/MD_Test_Dir
  ./install.sh --mode docker --db-profile test --root /data/MD_Test_Dir \
    --image metadetector:test --threads 8 --memory 16 --yes
