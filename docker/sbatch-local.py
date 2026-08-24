#!/usr/bin/env python3
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

VERSION = "MetaDetector local sbatch 0.1"

def directive(script, key):
    pat = re.compile(r"^\s*#SBATCH\s+" + re.escape(key) + r"(?:=|\s+)(.+?)\s*$")
    with open(script, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                return m.group(1).strip()
    return None

def all_directives(script):
    vals = {}
    with open(script, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = re.match(r"^\s*#SBATCH\s+(--[^=\s]+)(?:=|\s+)(.+?)\s*$", line)
            if m:
                vals[m.group(1)] = m.group(2).strip()
    return vals

if "--version" in sys.argv[1:]:
    print(VERSION)
    raise SystemExit(0)

args = sys.argv[1:]
script = None
for a in reversed(args):
    if not a.startswith("-"):
        script = a
        break

if not script:
    print("sbatch-local: no job script supplied", file=sys.stderr)
    raise SystemExit(2)

script = os.path.abspath(script)
if not os.path.isfile(script):
    print(f"sbatch-local: missing job script: {script}", file=sys.stderr)
    raise SystemExit(2)

state_dir = Path(os.environ.get("MD_LOCAL_SBATCH_STATE_DIR", "/tmp/md-sbatch-state"))
state_dir.mkdir(parents=True, exist_ok=True)
counter = state_dir / "next_job_id"
try:
    jobid = int(counter.read_text().strip())
except Exception:
    jobid = 1000
counter.write_text(str(jobid + 1))

directives = all_directives(script)
stdout_path = directives.get("--output")
stderr_path = directives.get("--error")
open_mode = directives.get("--open-mode", "truncate").lower()
dependency = directives.get("--dependency")

# The pipeline submits dependency metadata inside the job script.
# Because this shim executes jobs synchronously, parent jobs have already
# completed by the time a dependent job is submitted. Refuse execution if
# any recorded afterok dependency failed.
if dependency and dependency.startswith("afterok:"):
    deps = [x for x in dependency.split(":")[1:] if x]
    failed = []
    for dep in deps:
        rc_file = state_dir / f"{dep}.rc"
        if rc_file.exists():
            try:
                if int(rc_file.read_text().strip()) != 0:
                    failed.append(dep)
            except Exception:
                failed.append(dep)
        else:
            failed.append(dep)
    if failed:
        (state_dir / f"{jobid}.rc").write_text("1\n")
        print(f"Submitted batch job {jobid}")
        raise SystemExit(0)

mode = "ab" if open_mode == "append" else "wb"

def open_target(path):
    if not path:
        return subprocess.DEVNULL
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    return open(p, mode)

out_fh = open_target(stdout_path)
err_fh = open_target(stderr_path)

env = os.environ.copy()
env["SLURM_JOB_ID"] = str(jobid)
env["SLURM_JOB_NAME"] = directives.get("--job-name", f"md-{jobid}")
env["SLURM_CPUS_PER_TASK"] = directives.get("--cpus-per-task", "1")
env["SLURM_SUBMIT_DIR"] = os.getcwd()

try:
    proc = subprocess.run(
        ["/bin/bash", script],
        cwd=os.getcwd(),
        env=env,
        stdout=out_fh,
        stderr=err_fh,
        check=False,
    )
    rc = proc.returncode
finally:
    if hasattr(out_fh, "close"):
        out_fh.close()
    if hasattr(err_fh, "close"):
        err_fh.close()

(state_dir / f"{jobid}.rc").write_text(f"{rc}\n")
(state_dir / f"{jobid}.script").write_text(script + "\n")

# Slurm reports successful submission, not the eventual job exit status.
# Preserve that behavior so proc_assembly.sh receives a job ID.
print(f"Submitted batch job {jobid}")
raise SystemExit(0)
