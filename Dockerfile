FROM almalinux:9

SHELL ["/bin/bash", "-lc"]

RUN dnf install -y \
      bzip2 \
      ca-certificates \
      curl-minimal \
      findutils \
      gzip \
      procps-ng \
      tar \
      unzip \
      wget \
      which \
      xz \
      git \
    && dnf clean all \
    && rm -rf /var/cache/dnf

ARG MEGAN_VERSION=7_1_1

RUN mkdir -p /opt/megan /tmp/megan-install && \
    wget -O /tmp/megan-install/MEGAN_Community_unix_${MEGAN_VERSION}.sh \
      https://software-ab.cs.uni-tuebingen.de/download/megan7/MEGAN_Community_unix_${MEGAN_VERSION}.sh && \
    bash /tmp/megan-install/MEGAN_Community_unix_${MEGAN_VERSION}.sh \
      -q \
      -dir /opt/megan && \
    test -x /opt/megan/tools/daa-meganizer && \
    rm -rf /tmp/megan-install
    
ARG MAMBA_ROOT_PREFIX=/opt/micromamba
ENV MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX}

RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvj -C /tmp bin/micromamba \
    && install -m 0755 /tmp/bin/micromamba /usr/local/bin/micromamba \
    && rm -rf /tmp/bin


WORKDIR /opt

COPY md.yml /tmp/md.yml

RUN python3 - <<'PY_MD_YML'
from pathlib import Path

path = Path("/tmp/md.yml")
lines = path.read_text().splitlines()

out = []
in_channels = False
saw_channels = False
nodefaults_seen = False

for line in lines:
    stripped = line.strip()

    if stripped == "channels:":
        saw_channels = True
        in_channels = True
        out.append(line)
        continue

    if in_channels:
        if line and not line.startswith((" ", "\t")):
            if not nodefaults_seen:
                out.append("  - nodefaults")
            in_channels = False

        if in_channels:
            if stripped in ("- defaults", "defaults"):
                continue
            if stripped == "- nodefaults":
                nodefaults_seen = True

    out.append(line)

if in_channels and not nodefaults_seen:
    out.append("  - nodefaults")

if not saw_channels:
    out = [
        "channels:",
        "  - conda-forge",
        "  - bioconda",
        "  - nodefaults",
        *out,
    ]

path.write_text("\n".join(out) + "\n")
PY_MD_YML

RUN echo "===== Effective md.yml =====" \
    && cat /tmp/md.yml \
    && echo "============================" \
    && micromamba create -y \
         -n md \
         -f /tmp/md.yml \
         -c conda-forge \
         -c bioconda \
    && micromamba install -y \
         -n md \
         -c conda-forge \
         -c bioconda \
         mmseqs2 \
         pigz \
         blast \
         sra-tools \
    && micromamba clean --all --yes \
    && rm -f /tmp/md.yml

ENV PATH="${MAMBA_ROOT_PREFIX}/envs/md/bin:${PATH}"


COPY . /opt/metadetector

COPY docker/sbatch-local.py /usr/local/bin/sbatch
COPY docker/test-sbatch-local.sh /usr/local/bin/test-sbatch-local

RUN chmod 0755 \
      /usr/local/bin/sbatch \
      /usr/local/bin/test-sbatch-local \
      /opt/metadetector/queue.sh \
      /opt/metadetector/proc_assembly.sh \
    && find /opt/metadetector/helper_scripts /opt/metadetector/update_scripts \
         -type f \( -name '*.sh' -o -name '*.py' -o -name '*.pl' \) \
         -exec chmod a+rx {} + 2>/dev/null || true


RUN cat > /usr/local/bin/conda <<'EOF_CONDA_WRAPPER'
#!/usr/bin/env bash
set -Eeuo pipefail

case "${1:-}" in
    activate)
        env_name="${2:-}"
        if [[ "$env_name" != "md" ]]; then
            echo "ERROR: only the md environment is installed" >&2
            exit 1
        fi
        exit 0
        ;;
    run)
        shift
        exec micromamba run "$@"
        ;;
    env)
        shift
        exec micromamba env "$@"
        ;;
    list)
        exec micromamba list -n md "${@:2}"
        ;;
    *)
        exec micromamba "$@"
        ;;
esac
EOF_CONDA_WRAPPER

RUN chmod 0755 /usr/local/bin/conda

RUN mkdir -p \
      /database \
      /input \
      /output \
      /scratch \
      /state/sbatch \
    && chmod 0777 \
      /output \
      /scratch \
      /state \
      /state/sbatch

ENV MD_LOCAL_SBATCH_STATE_DIR=/state/sbatch
ENV TMPDIR=/scratch

WORKDIR /opt/metadetector


RUN python --version \
    && diamond version \
    && mmseqs version \
    && pigz --version \
    && makeblastdb -version \
    && fasterq-dump --version \
    && spades.py --version \
    && sbatch --version \
    && test-sbatch-local

CMD ["/bin/bash"]
