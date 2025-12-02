#!/bin/bash
# Pipeline to process metagenomic samples

source md.config
partition=${md_partition} # SLURM queue to process jobs
dbdir=${md_dbdir} # database directory
meganpath=${md_meganpath} # path to MEGAN tools
diamond_dbdir=${md_diamond_dbdir} # diamond database directory
diamond_args=${md_diamond_args} # memory args for diamond blastx
mmseqs_dbdir=${md_mmseqs_dbdir} # mmseqs nt database directory
mmseqs_tmpdir=${md_mmseqs_tmpdir} # mmseqs nt database directory
blastdb_nt=${dbdir}/blastdb/nt/nt
megandb=${dbdir}/megan
host_db_human=${md_human_ref}
host_db_contam=${md_contam_ref}
host_db_silva=${md_silva_ref}

#host_db_human="${dbdir}/bbmap/human_CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
#host_db_contam="${dbdir}/bbmap/silva_release132_oral_flora_rRNA/silva_release132_oral_flora_rRNA.fa"
#host_db_silva="${dbdir}/bbmap/silva_release132_oral_flora_rRNA/silva_release132_oral_flora_rRNA.fa"


eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate md
today=$(date +"%H:%M %d %b %Y")

display_usage() {
	echo -e "\nUsage:\n$0 [-a assembly_flag] [-j assembly_type] [-b BBMap args] [-d host_db for host removal] [-h help] [-k BBDuk args][-r path_to_sample] [-s stages] [-t threads] [-m memory]\n" >&2
}

while getopts "h?a:b:d:k:m:r:s:t:j:" opt; do
       	case "$opt" in
		a)	assembly_flag="${OPTARG}" ;;
	    b)	bbmap_args="${OPTARG}" ;;
	    d)	host_db="$OPTARG" ;;
	    h|\?)	display_usage
			exit 1
			;;
	    k)	bbduk_args="${OPTARG}" ;;
	    m)  memory="$OPTARG" ;;
	    r)  path="$OPTARG" ;;
	    s)  csv_stages="$OPTARG" ;;
	    t)  threads="$OPTARG" ;;
	    j)  assembly_type="$OPTARG" ;; # s (short illumina reads), h (hybrid with both short and long), l (long ONT reads)
       	esac
done

shift $((OPTIND-1))
sample=$(basename ${path})
sample_path=$(basename "$(dirname "${path}")")/$(basename "${path}")
account="MD:${sample_path}"


logpath=${path}/log
jobpath=${path}/job_script
statuspath=${path}/status_log
last_stage=26
spades_mem=${memory}
[[ ${spades_mem} -lt 250 ]] && spades_mem=250

assembly_flag=$(echo ${assembly_flag} | xargs) # remove whitespace
assembly_type=$(echo ${assembly_type} | xargs)
bbmap_args=$(echo ${bbmap_args} | xargs) # remove whitespace
bbduk_args=$(echo ${bbduk_args} | xargs) # remove whitespace
csv_stages=$(echo ${csv_stages} | xargs) # remove whitespace
host_db=$(echo ${host_db} | xargs) # remove whitespace
if [[ -z "${csv_stages}" ]]
then
       	[[ ${assembly_flag} == "a" && ${assembly_type} == "s" ]] && csv_stages=1,2,3a,4,5,6,7,9,10,12,13,14,16,17,18,19,21,22,23,25,26 # run all except metaviral assembly stages
       	[[ ${assembly_flag} == "i" && ${assembly_type} == "s" ]] && csv_stages=1,2,3a,4,5,7,10,13,14,16,18,19,21,23,25,26 # run for isolate only assembly
       	[[ ${assembly_flag} == "m" && ${assembly_type} == "s" ]] && csv_stages=1,2,3a,4,5,6,9,12,14,16,17,19,21,22,25,26 # run for metagenomic only assembly

	[[ ${assembly_flag} == "a" && ${assembly_type} == "l" ]] && csv_stages=1,2b,3c,4b,5b,6c,9c,12c,14b,17c,19b,23c,23cc,25,26 # run all except metaviral assembly stages
       	[[ ${assembly_flag} == "i" && ${assembly_type} == "l" ]] && csv_stages=1,2b,3c,4b,5b,6c,9c,12c,14b,17c,19b,23cc,25,26 # run for isolate only assembly
       	[[ ${assembly_flag} == "m" && ${assembly_type} == "l" ]] && csv_stages=1,2b,3c,4b,5b,6c,9c,12c,14b,17c,19b,23c,25,26 # run for metagenomic only assembly

	[[ ${assembly_flag} == "a" && ${assembly_type} == "h" ]] && csv_stages=1,2,2b,3a,3c,4,4b,5,5b,6b,7b,9b,9d,12b,12d,14,14b,17b,17d,19,19b,23b,23bb,23d,25,26 # run all except metaviral assembly stages
       	[[ ${assembly_flag} == "i" && ${assembly_type} == "h" ]] && csv_stages=1,2,2b,3a,3c,4,4b,5,5b,6b,7b,9b,9d,12b,12d,14,14b,17b,17d,19,19b,23bb,23d,25,26 # run for isolate only assembly
       	[[ ${assembly_flag} == "m" && ${assembly_type} == "h" ]] && csv_stages=1,2,2b,3a,3c,4,4b,5,5b,6b,7b,9b,9d,12b,12d,14,14b,17b,17d,19,19b,23b,25,26 # run for metagenomic only assembly

	[[ ${assembly_flag} == "t" && ${assembly_type} == "s" ]] && csv_stages=1,2,3b,4,5,6,7,9,10,12,13,14,16,17,18,19,21,22,23,25,26 # run all except metaviral assembly stages with targeted mapping filter
       	[[ ${assembly_flag} == "v" && ${assembly_type} == "s" ]] && csv_stages=1,2,3a,4,5,8,11,14,15,19,20,24,25,26 # run for metaviral only assembly
fi


[[ -z ${csv_stages} ]] && stages=( $(seq ${last_stage}) ) || IFS="," read -r -a stages <<< ${csv_stages}; unset IFS

echo "Assembly flag: $assembly_flag"
echo "Assembly type: $assembly_type"
echo "The stages for these parameters are: ${stages[@]}"

containsElement () {
	local e match="$1"
	shift
	for e; do [[ "$e" == "$match" ]] && return 0; done
	return 1
}


# set bbduk/bbmap args
[[ -z "${bbmap_args}" ]] && bbmap_args="fast usejni tossbrokenreads"
[[ -z "${bbduk_args}" ]] && bbduk_args="ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=rl trimq=30 maq=20 minlen=50 entropy=0.2 trimpolyg=10 trimpolya=10 ref=phix,adapters"
#echo bbduk_args="${bbduk_args}"

# set bbmap default script for host removal
bbmap="bbmap.sh"

# set host_db
[[ ! -z ${host_db} ]] && host_db="ref=${host_db} nodisk"
[[ -z ${host_db} ]] && host_db="ref=${host_db_human} nodisk"
[[ ! -z ${host_db_silva} ]] && host_db_silva="ref=${host_db_silva} nodisk"
[[ ! -z ${host_db_contam} ]] && host_db_contam="ref=${host_db_contam} nodisk"

echo " Human database: ${host_db_human}"



echo "*****************************"
echo "Starting process from STAGE ${stages[0]} at $today"
[ -d ${logpath} ] && echo "${logpath} exists. ! creating." \
       	|| mkdir -p ${logpath}
[ -d ${jobpath} ] && echo "${jobpath} exists. ! creating." \
       	|| mkdir -p ${jobpath}
[ -d ${statuspath} ] && echo "${statuspath} exists. ! creating." \
       	|| mkdir -p ${statuspath}
jobscript=()

# STAGE 1: produce FastQC pretrim reports
if containsElement "1" "${stages[@]}";
then
       	echo "STAGE 1: producing FastQC pretrim reports"
       	[ -d ${path}/fastqc/pretrim ] && \
		echo "${path}/fastqc/pretrim/shotgun exists. ! creating." \
		|| mkdir -p ${path}/fastqc/pretrim
	job=${jobpath}/stage1.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage1" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=2" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage1.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage1.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage1.finished ]]; then" >> ${job}
	echo -e "	if [[ ${assembly_type} == 's' ]] || [[ ${assembly_type} == 'h' ]]; then" >> ${job}
	echo -e "		fastqc --outdir ${path}/fastqc/pretrim --java java --threads=2 ${path}/${sample}_*R[12]*.fastq.gz" >> ${job}
	echo -e "       fi" >> ${job}
	echo -e "	if [[ ${assembly_type} == 'l' ]] || [[ ${assembly_type} == 'h' ]]; then" >> ${job}
	echo -e "		fastqc --outdir ${path}/fastqc/pretrim --java java --threads=1 ${path}/${sample}.fastq.gz" >> ${job}
        echo -e "       fi" >> ${job}
	echo -e "	out=\$?" >> ${job}
	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
	echo -e "		echo error \${out} running fastqc" >> ${job}
	echo -e "		touch ${statuspath}/stage1.not_finished" >> ${job}
	echo -e "	else" >> ${job}
	echo -e "		touch ${statuspath}/stage1.finished" >> ${job}
	echo -e "	fi" >> ${job}
	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage1=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage1}")
fi

# STAGE 2: quality trim shotgun libraries
if containsElement "2" "${stages[@]}";
then
       	echo "STAGE 2: quality trimming shotgun reads"
       	[ -d ${path}/trim/shotgun ] && \
		echo "${path}/trim/shotgun exists. ! creating." \
	       	|| mkdir -p ${path}/trim/shotgun
	job=${jobpath}/stage2.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage2" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage2.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage2.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "r1=${path}/${sample}_R1.fastq.gz" >> ${job}
	echo -e "r2=${path}/${sample}_R2.fastq.gz" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage2.finished ]]; then" >> ${job}
       	echo -e "	bbduk.sh "${bbduk_args}" ow \
			in1=\$r1 in2=\$r2 \
			out=${path}/trim/shotgun/${sample}_trimmed.fastq.gz \
			stats=${path}/trim/${sample}_bbduk_stats.txt \
			refstats=${path}/trim/${sample}_bbduk_refstats.txt" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running bbduk" >> ${job}
       	echo -e "		touch ${statuspath}/stage2.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
        echo -e "       	seqtk seq -1 ${path}/trim/shotgun/${sample}_trimmed.fastq.gz | \
	                	pigz -1 > ${path}/trim/shotgun/${sample}_trimmed_R1.fastq.gz; \
        	         	seqtk seq -2 ${path}/trim/shotgun/${sample}_trimmed.fastq.gz | \
                		pigz -1 > ${path}/trim/shotgun/${sample}_trimmed_R2.fastq.gz" >> ${job}
       	echo -e "		out=\$?" >> ${job}
       	echo -e	"		if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "			echo error \${out} running seqtk" >> ${job}
       	echo -e "			touch ${statuspath}/stage2.not_finished" >> ${job}
       	echo -e "		else" >> ${job}
       	echo -e "			touch ${statuspath}/stage2.finished" >> ${job}
       	echo -e "		fi" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage2=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage2}")
fi

# STAGE 2b: quality trim shotgun libraries
if containsElement "2b" "${stages[@]}";
then
       	echo "STAGE 2b: quality trimming ONT reads"
       	[ -d ${path}/trim/shotgun ] && \
		echo "${path}/trim/shotgun exists. ! creating." \
	       	|| mkdir -p ${path}/trim/shotgun
	job=${jobpath}/stage2b.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage2b" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage2b.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage2b.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "r=${path}/${sample}.fastq.gz" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage2b.finished ]]; then" >> ${job}
	echo -e "	fastp --dedup --overrepresentation_analysis --low_complexity_filter --trim_poly_x --report_title "${sample}_fastp_report" -q 10 -e 10 --cut_front --cut_tail -w 16 --in1 \${r} -o ${path}/trim/shotgun/${sample}_LR_trimmed.fastq.gz" >> ${job}
	echo -e "	out=\$?" >> ${job}
	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
	echo -e "		echo error \${out} running fastp" >> ${job}
	echo -e "		touch ${statuspath}/stage2b.not_finished" >> ${job}
	echo -e "	else" >> ${job}
	echo -e "		touch ${statuspath}/stage2b.finished" >> ${job}
	echo -e "	fi" >> ${job}
	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage2b=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage2b}")
fi

# STAGE 3a: remove provided host sequence from reads (dependency STAGE 2)
if containsElement "3a" "${stages[@]}";
then
       	echo "STAGE 3a: removing host from shotgun reads"
       	[ -d ${path}/trim/shotgun ] && \
		echo "${path}/trim/shotgun exists. ! creating." \
	       	|| mkdir -p ${path}/trim/shotgun
	job=${jobpath}/stage3.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage3a" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage3a.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage3a.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage2} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage2}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage3a.finished ]]; then" >> ${job}
       	echo -e "	${bbmap} in=${path}/trim/shotgun/${sample}_trimmed.fastq.gz \
			"${host_db}" "${bbmap_args}" covstats=${path}/trim/${sample}_host_covstats.txt \
			statsfile=${path}/trim/${sample}_host_statsfile.txt \
			ow 32bit pigz -Xmx${memory}g \
			outm=${path}/trim/shotgun/${sample}_host.fastq.gz \
			outu=${path}/trim/shotgun/${sample}_host_removed_pe.fastq.gz" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running bbmap.sh" >> ${job}
       	echo -e "		touch ${statuspath}/stage3a.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		seqtk seq -1 ${path}/trim/shotgun/${sample}_host_removed_pe.fastq.gz | \
	       			pigz -1 > ${path}/trim/shotgun/${sample}_host_removed_R1.fastq.gz; \
		 		seqtk seq -2 ${path}/trim/shotgun/${sample}_host_removed_pe.fastq.gz | \
				pigz -1 > ${path}/trim/shotgun/${sample}_host_removed_R2.fastq.gz" >> ${job}
       	echo -e "		out=\$?" >> ${job}
       	echo -e	"		if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "			echo error \${out} running seqtk" >> ${job}
       	echo -e "			touch ${statuspath}/stage3a.not_finished" >> ${job}
       	echo -e "		else" >> ${job}
       	echo -e "			touch ${statuspath}/stage3a.finished" >> ${job}
       	echo -e "		fi" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage3=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage3}")
fi

# STAGE 3b: map reads to reference sequence and use for downstream analysis (dependency STAGE 2)
if containsElement "3b" "${stages[@]}";
then
       	echo "STAGE 3b: keeping host from shotgun reads to use for downstream analysis"
       	[ -d ${path}/trim/shotgun ] && \
		echo "${path}/trim/shotgun exists. ! creating." \
	       	|| mkdir -p ${path}/trim/shotgun
	job=${jobpath}/stage3.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage3b" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage3b.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage3b.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage2} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage2}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage3b.finished ]]; then" >> ${job}
       	echo -e "	bbmap.sh in=${path}/trim/shotgun/${sample}_trimmed.fastq.gz \
			"${host_db}" "${bbmap_args}" covstats=${path}/trim/${sample}_host_covstats.txt \
			statsfile=${path}/trim/${sample}_host_statsfile.txt \
			ow 32bit pigz ambig=random -Xmx${memory}g \
			outu=${path}/trim/shotgun/${sample}_host.fastq.gz \
			outm=${path}/trim/shotgun/${sample}_host_removed_pe.fastq.gz" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running bbmap.sh" >> ${job}
       	echo -e "		touch ${statuspath}/stage3b.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		seqtk seq -1 ${path}/trim/shotgun/${sample}_host_removed_pe.fastq.gz | \
	       			pigz -1 > ${path}/trim/shotgun/${sample}_host_removed_R1.fastq.gz; \
		 		seqtk seq -2 ${path}/trim/shotgun/${sample}_host_removed_pe.fastq.gz | \
				pigz -1 > ${path}/trim/shotgun/${sample}_host_removed_R2.fastq.gz" >> ${job}
       	echo -e "		out=\$?" >> ${job}
       	echo -e	"		if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "			echo error \${out} running seqtk" >> ${job}
       	echo -e "			touch ${statuspath}/stage3b.not_finished" >> ${job}
       	echo -e "		else" >> ${job}
       	echo -e "			touch ${statuspath}/stage3b.finished" >> ${job}
       	echo -e "		fi" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage3=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage3}")
fi

# STAGE 3c: remove provided host sequence from reads (dependency STAGE 2)
if containsElement "3c" "${stages[@]}";
then
       	echo "STAGE 3c: removing host from ONT reads"
       	[ -d ${path}/trim/shotgun ] && \
		echo "${path}/trim/shotgun exists. ! creating." \
	       	|| mkdir -p ${path}/trim/shotgun
	job=${jobpath}/stage3c.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage3c" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage3c.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage3c.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage2b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage2b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage3c.finished ]]; then" >> ${job}
       	echo -e "	minimap2 -ax map-ont -t ${threads} "${host_db}" ${path}/trim/shotgun/${sample}_LR_trimmed.fastq.gz > ${path}/trim/shotgun/${sample}_host_removed_LR.sam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running minimap2" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		samtools fastq -f 4 ${path}/trim/shotgun/${sample}_host_removed_LR.sam > ${path}/trim/shotgun/${sample}_host_removed_LR.fastq" >> ${job}
       	echo -e "		samtools fastq -F 4 ${path}/trim/shotgun/${sample}_host_removed_LR.sam > ${path}/trim/shotgun/${sample}_host_LR.fastq" >> ${job}
       	echo -e "		touch ${statuspath}/stage3c.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage3c=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage3c}")
fi

# STAGE 4: remove common lab contaminant reads (dependency STAGE 3)
if containsElement "4" "${stages[@]}";
then
       	echo "STAGE 4: removing common lab contaminant reads"
       	[ -d ${path}/trim/shotgun ] && \
		echo "${path}/trim/shotgun exists. ! creating." \
	       	|| mkdir -p ${path}/trim/shotgun
	job=${jobpath}/stage4.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage4" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage4.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage4.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage3} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage3}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage4.finished ]]; then" >> ${job}
       	echo -e "	bbmap.sh in=${path}/trim/shotgun/${sample}_host_removed_pe.fastq.gz \
			"${host_db_contam}" tossbrokenreads printunmappedcount=t covstats=${path}/trim/${sample}_contaminant_covstats.txt \
			statsfile=${path}/trim/${sample}_contaminant_statsfile.txt \
			semiperfectmode=t maxindel=3 kfilter=25 maxsites=1 k=14 32bit pigz usejni -Xmx${memory}g \
			outm=${path}/trim/shotgun/${sample}_contaminant.fastq.gz \
			outu=${path}/trim/shotgun/${sample}_host_contaminant_removed_pe.fastq.gz overwrite=true" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running bbmap.sh" >> ${job}
       	echo -e "		touch ${statuspath}/stage4.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		seqtk seq -1 ${path}/trim/shotgun/${sample}_host_contaminant_removed_pe.fastq.gz | \
				pigz -1 > ${path}/trim/shotgun/${sample}_host_contaminant_removed_R1.fastq.gz; \
				seqtk seq -2 ${path}/trim/shotgun/${sample}_host_contaminant_removed_pe.fastq.gz | \
				pigz -1 > ${path}/trim/shotgun/${sample}_host_contaminant_removed_R2.fastq.gz" >> ${job}
       	echo -e "		out=\$?" >> ${job}
       	echo -e	"		if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "			echo error \${out} running seqtk" >> ${job}
       	echo -e "			touch ${statuspath}/stage4.not_finished" >> ${job}
       	echo -e "		else" >> ${job}
       	echo -e "			touch ${statuspath}/stage4.finished" >> ${job}
       	echo -e "		fi" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage4=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage4}")
fi

# STAGE 4b: remove common lab contaminant reads (dependency STAGE 3)
if containsElement "4b" "${stages[@]}";
then
       	echo "STAGE 4b: removing common lab contaminant reads from LR"
       	[ -d ${path}/trim/shotgun ] && \
		echo "${path}/trim/shotgun exists. ! creating." \
	       	|| mkdir -p ${path}/trim/shotgun
	job=${jobpath}/stage4b.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage4b" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage4b.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage4b.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage3c} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage3c}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage4b.finished ]]; then" >> ${job}
       	echo -e "	minimap2 -ax map-ont -N 1 -k 14 -t ${threads} "${host_db_contam}" ${path}/trim/shotgun/${sample}_host_removed_LR.fastq > ${path}/trim/shotgun/${sample}_host_contaminant_removed_LR.sam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running minimap2" >> ${job}
       	echo -e	"	else" >> ${job}
       	echo -e "		samtools fastq -f 4 ${path}/trim/shotgun/${sample}_host_contaminant_removed_LR.sam > ${path}/trim/shotgun/${sample}_host_contaminant_removed_LR.fastq" >> ${job}
       	echo -e "		samtools fastq -F 4 ${path}/trim/shotgun/${sample}_host_contaminant_removed_LR.sam > ${path}/trim/shotgun/${sample}_contaminant_LR.fastq" >> ${job}
       	echo -e "		touch ${statuspath}/stage4b.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage4b=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage4b}")
fi

# STAGE 5: remove common flora rRNA reads (dependency STAGE 4)
if containsElement "5" "${stages[@]}";
then
       	echo "STAGE 5: removing common flora rRNA reads"
       	[ -d ${path}/trim/shotgun ] && \
		echo "${path}/trim/shotgun exists. ! creating." \
	       	|| mkdir -p ${path}/trim/shotgun
	job=${jobpath}/stage5.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage5" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage5.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage5.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage4} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage4}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage5.finished ]]; then" >> ${job}
       	echo -e "	bbmap.sh "${bbmap_args}" in=${path}/trim/shotgun/${sample}_host_contaminant_removed_pe.fastq.gz \
			"${host_db_silva}" tossbrokenreads printunmappedcount=t covstats=${path}/trim/rRNA_covstats.txt \
			semiperfectmode=t \
			outm=${path}/trim/shotgun/${sample}_rRNA.fastq.gz usejni=t -Xmx${memory}g \
			outu=${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz overwrite=true" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running remove_contaminants.sh" >> ${job}
       	echo -e "		touch ${statuspath}/stage5.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		seqtk seq -1 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz | \
				pigz -1 > ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_R1.fastq.gz; \
				seqtk seq -2 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz | \
				pigz -1 > ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_R2.fastq.gz" >> ${job}
       	echo -e "		out=\$?" >> ${job}
       	echo -e	"		if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "			echo error \${out} running seqtk" >> ${job}
       	echo -e "			touch ${statuspath}/stage5.not_finished" >> ${job}
       	echo -e "		else" >> ${job}
       	echo -e "			touch ${statuspath}/stage5.finished" >> ${job}
       	echo -e "		fi" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage5=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage5}")
fi

# STAGE 5b: remove common flora rRNA reads (dependency STAGE 4)
if containsElement "5b" "${stages[@]}";
then
       	echo "STAGE 5b: removing common flora rRNA reads"
       	[ -d ${path}/trim/shotgun ] && \
		echo "${path}/trim/shotgun exists. ! creating." \
	       	|| mkdir -p ${path}/trim/shotgun
	job=${jobpath}/stage5b.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage5b" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage5b.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage5b.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage4b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage4b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage5b.finished ]]; then" >> ${job}
       	echo -e "	minimap2 -ax map-ont -t ${threads} "${host_db_silva}" ${path}/trim/shotgun/${sample}_host_contaminant_removed_LR.fastq > ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.sam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running removeflora.sh" >> ${job}
       	echo -e	"	else" >> ${job}
       	echo -e "		samtools fastq -f 4 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.sam > ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.fastq" >> ${job}
       	echo -e "		samtools fastq -F 4 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.sam > ${path}/trim/shotgun/${sample}_rRNA_LR.fastq" >> ${job}
       	echo -e "		touch ${statuspath}/stage5b.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage5b=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage5b}")
fi

# STAGE 6: metaSPAdes PE trimmed assembly (dependency STAGE 5)
if containsElement "6" "${stages[@]}";
then
       	echo "STAGE 6: metaSPAdes PE trimmed assembly"
       	[ -d ${path}/spades/meta_pe_trim ] && \
		echo "${path}/spades/meta_pe_trim exists. ! creating." \
		|| mkdir -p ${path}/spades/meta_pe_trim
	job=${jobpath}/stage6.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage6" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage6.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage6.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage5} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage5}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage6.finished ]]; then" >> ${job}
       	echo -e "	spades.py -t ${threads} --only-assembler --meta --memory ${spades_mem} \
			--pe1-12 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
			-o ${path}/spades/meta_pe_trim" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running metaSPAdes" >> ${job}
       	echo -e "		touch ${statuspath}/stage6.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage6.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage6=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage6}")
fi

#STAGE 6b: Dragonflye meta hybrid assembly (dependency STAGE 5)
if containsElement "6b" "${stages[@]}";
then
       	echo "STAGE 6b: dragonflye hybrid assembly"
       	[ -d ${path}/dragonflye/hybrid ] && \
		echo "${path}/dragonflye/meta exists. ! creating." \
		|| mkdir -p ${path}/dragonflye/hybrid
	job=${jobpath}/stage6b.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage6b" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage6b.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage6b.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage5b} && ${jobid_stage5} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage5b}:${jobid_stage5}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage6b.finished ]]; then" >> ${job}
	echo -e "	dragonflye --force --cpus ${threads} --tmpdir /dev/shm --ram ${memory} \
		--racon 4 --outdir ${path}/dragonflye/hybrid --reads ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.fastq \
		--R1 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_R1.fastq.gz --R2 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_R2.fastq.gz">> ${job}
	echo -e "	out=\$?" >> ${job}
	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
	echo -e "		echo error \${out} running dragonflye" >> ${job}
	echo -e "		touch ${statuspath}/stage6b.not_finished" >> ${job}
	echo -e "	else" >> ${job}
	echo -e "		touch ${statuspath}/stage6b.finished" >> ${job}
	echo -e "	fi" >> ${job}
	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage6b=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage6b}")
fi

#STAGE 6c: Dragonflye isolate long assembly (dependency STAGE 5)
if containsElement "6c" "${stages[@]}";
then
       	echo "STAGE 6c: dragonflye LR assembly"
       	[ -d ${path}/dragonflye/LR ] && \
		echo "${path}/dragonflye/LR exists. ! creating." \
		|| mkdir -p ${path}/dragonflye/LR
	job=${jobpath}/stage6c.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage6c" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage6c.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage6c.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage5b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage5b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage6c.finished ]]; then" >> ${job}
	echo -e "	dragonflye --force --cpus ${threads} --tmpdir /dev/shm --ram ${memory} \
			--racon 4 --outdir ${path}/dragonflye/LR --reads ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.fastq" >> ${job}
	echo -e "	out=\$?" >> ${job}
	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
	echo -e "		echo error \${out} running dragonflye" >> ${job}
	echo -e "		touch ${statuspath}/stage6c.not_finished" >> ${job}
	echo -e "	else" >> ${job}
	echo -e "		touch ${statuspath}/stage6c.finished" >> ${job}
	echo -e "	fi" >> ${job}
	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage6c=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage6c}")
fi

# STAGE 7: SPAdes PE trimmed assembly (dependency STAGE 5)
if containsElement "7" "${stages[@]}";
then
        echo "STAGE 7: SPAdes PE trimmed assembly"
        [ -d ${path}/spades/pe_trim ] && \
                echo "${path}/spades/pe_trim exists. ! creating." \
                || mkdir -p ${path}/spades/pe_trim
        job=${jobpath}/stage7.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage7" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage7.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage7.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage5} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage5}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage7.finished ]]; then" >> ${job}
        echo -e "       spades.py --isolate \
			--threads ${threads} --memory ${spades_mem} \
                        --pe1-12 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
                        -o ${path}/spades/pe_trim" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running SPAdes" >> ${job}
        echo -e "               touch ${statuspath}/stage7.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage7.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage7=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage7}")
fi

# STAGE 7b: unicycler isolate hybrid assembly (dependency STAGE 5)
if containsElement "7b" "${stages[@]}";
then
        echo "STAGE 7b: unicycler isolate hybrid assembly"
        [ -d ${path}/unicycler/isolate ] && \
                echo "${path}/unicycler/isolate exists. ! creating." \
                || mkdir -p ${path}/unicycler/isolate
        job=${jobpath}/stage7b.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage7b" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage7b.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage7b.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage5} && ${jobid_stage5b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage5}:${jobid_stage5b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage7b.finished ]]; then" >> ${job}
        echo -e "	unicycler -1 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_R1.fastq.gz \
        -2 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_R2.fastq.gz \
        -l ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.fastq --min_fasta_length 1000 \
		--threads ${threads} --mode normal --out ${path}/unicycler/isolate" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running SPAdes" >> ${job}
        echo -e "               touch ${statuspath}/stage7b.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage7b.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage7b=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage7b}")
fi

# STAGE 8: metaviralSPAdes PE trimmed assembly (dependency STAGE 5)
if containsElement "8" "${stages[@]}";
then
        echo "STAGE 8: metaviralSPAdes PE trimmed assembly"
        [ -d ${path}/spades/metaviral_pe_trim ] && \
                echo "${path}/spades/metaviral_pe_trim exists. ! creating." \
                || mkdir -p ${path}/spades/metaviral_pe_trim
        job=${jobpath}/stage8.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage8" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage8.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage8.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage5} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage5}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage8.finished ]]; then" >> ${job}
       	echo -e "	[[ ${memory} -lt 250 ]] && ${memory} = 250" >> ${job}
        echo -e "       spades.py --metaviral -t ${threads} --memory ${memory} \
                        --pe1-12 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
                        -o ${path}/spades/metaviral_pe_trim" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running SPAdes" >> ${job}
        echo -e "               touch ${statuspath}/stage8.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage8.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage8=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage8}")
fi

# STAGE 9: Map trim reads to trim metaSPAdes contig assembly (dependency STAGE 6)
if containsElement "9" "${stages[@]}";
then
       	echo "STAGE 9: Map trim reads to trim metaSPAdes contigs"
       	[ -d ${path}/mapping/meta_pe_trim ] && \
		echo "${path}/mapping/trim \
		exists. ! creating." || mkdir -p \
		${path}/mapping/meta_pe_trim
	job=${jobpath}/stage9.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage9" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage9.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage9.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage6} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage6}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage9.finished ]]; then" >> ${job}
       	echo -e "	bowtie2-build -q --threads ${threads} \
			-f ${path}/spades/meta_pe_trim/contigs.fasta ${path}/mapping/meta_pe_trim/${sample}" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running bowtie2" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	bowtie2 --met-file ${path}/mapping/meta_pe_trim/${sample}.met \
			--threads ${threads} -x ${path}/mapping/meta_pe_trim/${sample} \
			--interleaved ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
			--sensitive -S ${path}/mapping/meta_pe_trim/${sample}.sam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running bowtie2" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	samtools faidx ${path}/spades/meta_pe_trim/contigs.fasta" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running samtools" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	samtools view -@ ${threads} \
			-bt ${path}/spades/meta_pe_trim/contigs.fasta.fai ${path}/mapping/meta_pe_trim/${sample}.sam | \
			samtools sort -@ ${threads} - \
			-o ${path}/mapping/meta_pe_trim/${sample}_sorted.bam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running samtools" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	samtools index -@ ${threads} ${path}/mapping/meta_pe_trim/${sample}_sorted.bam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running samtools" >> ${job}
       	echo -e "		touch ${statuspath}/stage9.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		rm ${path}/mapping/meta_pe_trim/${sample}.sam" >> ${job}
       	echo -e "		touch ${statuspath}/stage9.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage9=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage9}")
fi

# STAGE 9b: Map trim reads to trim dragonflye contig assembly (dependency STAGE 6)
if containsElement "9b" "${stages[@]}";
then
       	echo "STAGE 9b: Map trim reads to trim dragonflye hybrid contigs"
       	[ -d ${path}/mapping/dragonflye_hyb_LR_trim ] && \
		echo "${path}/mapping/trim \
		exists. ! creating." || mkdir -p \
		${path}/mapping/dragonflye_hyb_LR_trim
	job=${jobpath}/stage9b.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage9b" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage9b.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage9b.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	echo -e "#SBATCH --dependency=afterok:${jobid_stage6b}" >> ${job}

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage9b.finished ]]; then" >> ${job}
       	echo -e "	minimap2 -ax map-ont ${path}/dragonflye/hybrid/contigs.fa ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.fastq > ${path}/mapping/dragonflye_hyb_LR_trim/${sample}.sam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running minimap2" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	samtools faidx ${path}/dragonflye/hybrid/contigs.fa" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running samtools" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	samtools view -@ ${threads} \
			-bt ${path}/dragonflye/hybrid/contigs.fa.fai ${path}/mapping/dragonflye_hyb_LR_trim/${sample}.sam | \
			samtools sort -@ ${threads} - \
			-o ${path}/mapping/dragonflye_hyb_LR_trim/${sample}_sorted.bam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running samtools" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	samtools index -@ ${threads} ${path}/mapping/dragonflye_hyb_LR_trim/${sample}_sorted.bam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running samtools" >> ${job}
       	echo -e "		touch ${statuspath}/stage9b.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		rm ${path}/mapping/dragonflye_hyb_LR_trim/${sample}.sam" >> ${job}
       	echo -e "		touch ${statuspath}/stage9b.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage9b=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage9b}")
fi

# STAGE 9c: Map trim reads to trim dragonflye contig assembly (dependency STAGE 6)
if containsElement "9c" "${stages[@]}";
then
       	echo "STAGE 9c: Map trim reads to trim dragonflye contigs"
       	[ -d ${path}/mapping/dragonflye_LR_trim ] && \
		echo "${path}/mapping/trim \
		exists. ! creating." || mkdir -p \
		${path}/mapping/dragonflye_LR_trim
	job=${jobpath}/stage9c.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage9c" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage9c.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage9c.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	echo -e "#SBATCH --dependency=afterok:${jobid_stage6c}" >> ${job}

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage9c.finished ]]; then" >> ${job}
       	echo -e "	minimap2 -ax map-ont ${path}/dragonflye/LR/contigs.fa ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.fastq > ${path}/mapping/dragonflye_LR_trim/${sample}.sam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running minimap2" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	samtools faidx ${path}/dragonflye/LR/contigs.fa" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running samtools" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	samtools view -@ ${threads} \
			-bt ${path}/dragonflye/LR/contigs.fa.fai ${path}/mapping/dragonflye_LR_trim/${sample}.sam | \
			samtools sort -@ ${threads} - \
			-o ${path}/mapping/dragonflye_LR_trim/${sample}_sorted.bam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running samtools" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	samtools index -@ ${threads} ${path}/mapping/dragonflye_LR_trim/${sample}_sorted.bam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running samtools" >> ${job}
       	echo -e "		touch ${statuspath}/stage9c.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		rm ${path}/mapping/dragonflye_LR_trim/${sample}.sam" >> ${job}
       	echo -e "		touch ${statuspath}/stage9c.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage9c=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage9c}")
fi

# STAGE 9d: Map trim reads to trim unicycler contig assembly (dependency STAGE 6)
if containsElement "9d" "${stages[@]}";
then
       	echo "STAGE 9d: Map trim reads to trim unicycler contigs"
       	[ -d ${path}/mapping/unicycler_hyb_trim ] && \
		echo "${path}/mapping/trim \
		exists. ! creating." || mkdir -p \
		${path}/mapping/unicycler_hyb_trim
	job=${jobpath}/stage9d.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage9d" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage9d.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage9d.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	echo -e "#SBATCH --dependency=afterok:${jobid_stage7b}" >> ${job}

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage9d.finished ]]; then" >> ${job}
       	echo -e "	minimap2 -ax map-ont ${path}/unicycler/isolate/assembly.fasta ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_isolate.fastq > ${path}/mapping/unicycler_hyb_trim/${sample}.sam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running minimap2" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	samtools faidx ${path}/unicycler/isolate/contigs.fasta" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running samtools" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	samtools view -@ ${threads} \
			-bt ${path}/unicycler/isolate/assembly.fasta.fai ${path}/mapping/unicycler_hyb_trim/${sample}.sam | \
			samtools sort -@ ${threads} - \
			-o ${path}/mapping/unicycler_hyb_trim/${sample}_sorted.bam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running samtools" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	samtools index -@ ${threads} ${path}/mapping/unicycler_hyb_trim/${sample}_sorted.bam" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running samtools" >> ${job}
       	echo -e "		touch ${statuspath}/stage9d.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		rm ${path}/mapping/unicycler_hyb_trim/${sample}.sam" >> ${job}
       	echo -e "		touch ${statuspath}/stage9d.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage9d=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage9d}")
fi

# STAGE 10: Map trim reads to trim SPAdes contig assembly (dependency STAGE 7)
if containsElement "10" "${stages[@]}";
then
        echo "STAGE 10: Map trim reads to trim SPAdes contigs"
        [ -d ${path}/mapping/pe_trim ] && \
                echo "${path}/mapping/pe_trim \
                exists. ! creating." || mkdir -p \
                ${path}/mapping/pe_trim
        job=${jobpath}/stage10.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage10" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage10.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage10.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage7} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage7}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage10.finished ]]; then" >> ${job}
        echo -e "       bowtie2-build -q --threads ${threads} \
                        -f ${path}/spades/pe_trim/contigs.fasta ${path}/mapping/pe_trim/${sample}" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running bowtie2" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "       bowtie2 --met-file ${path}/mapping/pe_trim/${sample}.met \
                        --threads ${threads} -x ${path}/mapping/pe_trim/${sample} \
                        --interleaved ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
                        --sensitive -S ${path}/mapping/pe_trim/${sample}.sam" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running bowtie2" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "       samtools faidx ${path}/spades/pe_trim/contigs.fasta" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running samtools" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "       samtools view -@ ${threads} \
                        -bt ${path}/spades/pe_trim/contigs.fasta.fai ${path}/mapping/pe_trim/${sample}.sam | \
                        samtools sort -@ ${threads} - \
                        -o ${path}/mapping/pe_trim/${sample}_sorted.bam" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running samtools" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "       samtools index -@ ${threads} ${path}/mapping/pe_trim/${sample}_sorted.bam" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running samtools" >> ${job}
        echo -e "               touch ${statuspath}/stage10.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               rm ${path}/mapping/pe_trim/${sample}.sam" >> ${job}
        echo -e "               touch ${statuspath}/stage10.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage10=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage10}")
fi

# STAGE 11: Map trim reads to trim metaviral SPAdes contig assembly (dependency STAGE 8)
if containsElement "11" "${stages[@]}";
then
        echo "STAGE 11: Map trim reads to trim SPAdes contigs"
        [ -d ${path}/mapping/metaviral_pe_trim ] && \
                echo "${path}/mapping/metaviral_pe_trim \
                exists. ! creating." || mkdir -p \
                ${path}/mapping/metaviral_pe_trim
        job=${jobpath}/stage11.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage11" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage11.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage11.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage8} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage8}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage11.finished ]]; then" >> ${job}
        echo -e "       bowtie2-build -q --threads ${threads} \
                        -f ${path}/spades/metaviral_pe_trim/contigs.fasta ${path}/mapping/metaviral_pe_trim/${sample}" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running bowtie2" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "       bowtie2 --met-file ${path}/mapping/metaviral_pe_trim/${sample}.met \
                        --threads ${threads} -x ${path}/mapping/metaviral_pe_trim/${sample} \
                        --interleaved ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
                        --sensitive -S ${path}/mapping/metaviral_pe_trim/${sample}.sam" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running bowtie2" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "       samtools faidx ${path}/spades/metaviral_pe_trim/contigs.fasta" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running samtools" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "       samtools view -@ ${threads} \
                        -bt ${path}/spades/metaviral_pe_trim/contigs.fasta.fai ${path}/mapping/metaviral_pe_trim/${sample}.sam | \
                        samtools sort -@ ${threads} - \
                        -o ${path}/mapping/metaviral_pe_trim/${sample}_sorted.bam" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running samtools" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "       samtools index -@ ${threads} ${path}/mapping/metaviral_pe_trim/${sample}_sorted.bam" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running samtools" >> ${job}
        echo -e "               touch ${statuspath}/stage11.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               rm ${path}/mapping/metaviral_pe_trim/${sample}.sam" >> ${job}
        echo -e "               touch ${statuspath}/stage11.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage11=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage11}")
fi

# STAGE 12: blastx metaSPAdes contigs (dependency STAGE 6)
if containsElement "12" "${stages[@]}";
then
       	echo "STAGE 12: diamond blastx metaSPAdes contigs"
       	[ -d ${path}/blast ] && echo "${path}/blast \
		exists. ! creating." || mkdir -p ${path}/blast
	job=${jobpath}/stage12.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage12" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage12.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage12.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage6} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage6}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage12.finished ]]; then" >> ${job}
       	echo -e "	diamond blastx ${diamond_args} --threads ${threads} \
			--db ${diamond_dbdir} --query ${path}/spades/meta_pe_trim/contigs.fasta \
			--range-culling -F 15 \
			--evalue 1e-5 \
			--outfmt 100 --out ${path}/blast/${sample}_metaspades_contigs_blastx.daa" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running diamond" >> ${job}
       	echo -e "		touch ${statuspath}/stage12.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage12.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage12=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage12}")
fi

# STAGE 12b: blastx dragonflye hybrid contigs (dependency STAGE 6)
if containsElement "12b" "${stages[@]}";
then
       	echo "STAGE 12b: diamond blastx dragonflye hybrid contigs"
       	[ -d ${path}/blast ] && echo "${path}/blast \
		exists. ! creating." || mkdir -p ${path}/blast
	job=${jobpath}/stage12b.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage12b" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage12b.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage12b.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage6b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage6b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage12b.finished ]]; then" >> ${job}
       	echo -e "	diamond blastx ${diamond_args} --threads ${threads} \
			--db ${diamond_dbdir} --query ${path}/dragonflye/hybrid/contigs.fa \
			--range-culling -F 15 \
			--evalue 1e-5 \
			--outfmt 100 --out ${path}/blast/${sample}_dragonflye_hybrid_contigs_blastx.daa" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running diamond" >> ${job}
       	echo -e "		touch ${statuspath}/stage12b.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage12b.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage12b=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage12b}")
fi

# STAGE 12c: blastx dragonflye LR contigs (dependency STAGE 6)
if containsElement "12c" "${stages[@]}";
then
       	echo "STAGE 12c: diamond blastx dragonflye LR contigs"
       	[ -d ${path}/blast ] && echo "${path}/blast \
		exists. ! creating." || mkdir -p ${path}/blast
	job=${jobpath}/stage12c.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage12c" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage12c.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage12c.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage6c} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage6c}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage12c.finished ]]; then" >> ${job}
       	echo -e "	diamond blastx ${diamond_args} --threads ${threads} \
			--db ${diamond_dbdir} --query ${path}/dragonflye/LR/contigs.fa \
			--range-culling -F 15 \
			--evalue 1e-5 \
			--outfmt 100 --out ${path}/blast/${sample}_dragonflye_LR_contigs_blastx.daa" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running diamond" >> ${job}
       	echo -e "		touch ${statuspath}/stage12c.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage12c.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage12c=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage12c}")
fi

# STAGE 12d: blastx unicycler contigs (dependency STAGE 6)
if containsElement "12d" "${stages[@]}";
then
       	echo "STAGE 12d: diamond blastx unicycler isolate contigs"
       	[ -d ${path}/blast ] && echo "${path}/blast \
		exists. ! creating." || mkdir -p ${path}/blast
	job=${jobpath}/stage12d.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage12d" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage12d.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage12d.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage7b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage7b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

	echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage12d.finished ]]; then" >> ${job}
       	echo -e "	diamond blastx ${diamond_args} --threads ${threads} \
			--db ${diamond_dbdir} --query ${path}/unicycler/isolate/assembly.fasta \
			--range-culling -F 15 \
			--evalue 1e-5 \
			--outfmt 100 --out ${path}/blast/${sample}_unicycler_contigs_blastx.daa" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running diamond" >> ${job}
       	echo -e "		touch ${statuspath}/stage12d.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage12d.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage12d=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage12d}")
fi

# STAGE 13: blastx SPAdes contigs (dependency STAGE 7)
if containsElement "13" "${stages[@]}";
then
        echo "STAGE 13: diamond blastx SPAdes contigs"
        [ -d ${path}/blast ] && echo "${path}/blast \
                exists. ! creating." || mkdir -p ${path}/blast
        job=${jobpath}/stage13.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage13" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage13.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage13.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage7} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage7}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage13.finished ]]; then" >> ${job}
       	echo -e "	diamond blastx ${diamond_args} --threads ${threads} \
                        --db ${diamond_dbdir} --query ${path}/spades/pe_trim/contigs.fasta \
			--range-culling -F 15 \
			--evalue 1e-5 \
                        --outfmt 100 --out ${path}/blast/${sample}_spades_contigs_blastx.daa" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running diamond" >> ${job}
        echo -e "               touch ${statuspath}/stage13.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage13.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage13=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage13}")
fi

# STAGE 14: blastx trimmed reads (dependency STAGE 5)
if containsElement "14" "${stages[@]}";
then
       	echo "STAGE 14: diamond blastx trimmed reads"
       	[ -d ${path}/blast ] && echo "${path}/blast \
		exists. ! creating." || mkdir -p ${path}/blast
	job=${jobpath}/stage14.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage14" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage14.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage14.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage5} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage5}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage14.finished ]]; then" >> ${job}
       	echo -e "	diamond blastx ${diamond_args} --threads ${threads} \
			--db ${diamond_dbdir} --query ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
			--evalue 1e-5 \
			--outfmt 100 --out ${path}/blast/${sample}_reads_blastx.daa" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running diamond" >> ${job}
       	echo -e "		touch ${statuspath}/stage14.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage14.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage14=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage14}")
fi

# STAGE 14b: blastx trimmed long reads (dependency STAGE 5)
if containsElement "14b" "${stages[@]}";
then
       	echo "STAGE 14b: diamond blastx trimmed reads"
       	[ -d ${path}/blast ] && echo "${path}/blast \
		exists. ! creating." || mkdir -p ${path}/blast
	job=${jobpath}/stage14b.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage14b" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage14b.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage14b.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage5b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage5b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage14b.finished ]]; then" >> ${job}
       	echo -e "	diamond blastx ${diamond_args} --threads ${threads} \
			--db ${diamond_dbdir} --query ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.fastq \
			--evalue 1e-5 \
			--outfmt 100 --out ${path}/blast/${sample}_LR_blastx.daa" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running diamond" >> ${job}
       	echo -e "		touch ${statuspath}/stage14b.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage14b.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage14b=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage14b}")
fi

# STAGE 15: blastx metaviral SPAdes contigs (dependency STAGE 8)
if containsElement "15" "${stages[@]}";
then
        echo "STAGE 15: diamond blastx metaviral SPAdes contigs"
        [ -d ${path}/blast ] && echo "${path}/blast \
                exists. ! creating." || mkdir -p ${path}/blast
        job=${jobpath}/stage15.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage15" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage15.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage15.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage8} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage8}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage15.finished ]]; then" >> ${job}
       	echo -e "	diamond blastx ${diamond_args} --threads ${threads} \
                        --db ${diamond_dbdir} --query ${path}/spades/metaviral_pe_trim/contigs.fasta \
			--evalue 1e-5 \
                        --outfmt 100 --out ${path}/blast/${sample}_metaviral_contigs_blastx.daa" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running diamond" >> ${job}
        echo -e "               touch ${statuspath}/stage15.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage15.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage15=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage15}")
fi

# STAGE 16: mmseqs2 metaSPAdes contigs (dependency STAGE 6)
if containsElement "16" "${stages[@]}";
then
       	echo "STAGE 16: mmseqs2 metaSPAdes contigs"
       	[ -d ${path}/blast ] && echo "${path}/blast \
		exists. ! creating." || mkdir -p ${path}/blast
	job=${jobpath}/stage16.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage16" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage16.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage16.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage6} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage6}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage16.finished ]]; then" >> ${job}
       	echo -e "	export TTY=0; mmseqs easy-search ${path}/spades/meta_pe_trim/contigs.fasta ${mmseqs_dbdir} \
		       	${path}/blast/${sample}_metaspades_contigs.mmseqs.out ${mmseqs_tmpdir} \
			--max-accept 10 --format-mode 0 --split-mode 0 \
			-s 5.7 -e 1.0e-5 --search-type 3 \
			--split-memory-limit ${memory}G --local-tmp ${mmseqs_tmpdir} --threads ${threads}" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running mmseqs easy-search" >> ${job}
       	echo -e "		touch ${statuspath}/stage16.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage16.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage16=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage16}")
fi

# STAGE 16: megablast metaSPAdes contigs (dependency STAGE 6)
if containsElement "X" "${stages[@]}";
then
       	echo "STAGE 16: blastn metaSPAdes contigs"
       	[ -d ${path}/blast ] && echo "${path}/blast \
		exists. ! creating." || mkdir -p ${path}/blast
	job=${jobpath}/stage16.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage16" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage16.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage16.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage6} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage6}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage16.finished ]]; then" >> ${job}
       	echo -e "	blastn -num_threads ${threads} -max_target_seqs 10 \
			-db ${blastdb_nt} -query ${path}/spades/meta_pe_trim/contigs.fasta \
			-evalue 1e-5 -outfmt 6 -out ${path}/blast/${sample}_metaspades_contigs.blastn" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running megablast" >> ${job}
       	echo -e "		touch ${statuspath}/stage16.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage16.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage16=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage16}")
fi

# STAGE 17: meganize blastx metaSPAdes contig output (dependency STAGE 12)
if containsElement "17" "${stages[@]}";
then
       	echo "STAGE 17: meganize blastx metaSPAdes contig output"
       	[ -d ${path}/blast ] && echo "${path}/blast \
		exists. ! creating." || mkdir -p ${path}/blast
	job=${jobpath}/stage17.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage17" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --chdir=${megandb}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage17.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage17.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage12} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage12}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage17.finished ]]; then" >> ${job}
	echo -e "	${meganpath}/daa-meganizer --in ${path}/blast/${sample}_metaspades_contigs_blastx.daa \
			--mapDB ${megandb}/megan-map.db --acc2taxa ${megandb}/ncbi.map --threads ${threads} --topPercent 0.5 \
			--minSupportPercent 0 --lcaAlgorithm longReads --longReads true --verbose" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running daa-meganizer" >> ${job}
       	echo -e "		touch ${statuspath}/stage17.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage17.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	paste <(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_metaspades_contigs_blastx.daa -c2c Taxonomy | awk '{print \$1}') \
				<(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_metaspades_contigs_blastx.daa -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') \
				> ${path}/blast/${sample}_metaspades_contigs_blastx_daa_summary_count.tsv" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage17=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage17}")
fi

# STAGE 17b: meganize blastx dragonflye contig output (dependency STAGE 15)
if containsElement "17b" "${stages[@]}";
then
        echo "STAGE 17b: meganize blastx dragonflye hybrid contig output"
        [ -d ${path}/blast ] && echo "${path}/blast \
                exists. ! creating." || mkdir -p ${path}/blast
        job=${jobpath}/stage17b.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage17b" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --chdir=${megandb}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage17b.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage17b.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage12b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage12b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage17b.finished ]]; then" >> ${job}
        echo -e "       ${meganpath}/daa-meganizer --in ${path}/blast/${sample}_dragonflye_hybrid_contigs_blastx.daa \
			--mapDB ${megandb}/megan-map.db --acc2taxa ${megandb}/ncbi.map --threads ${threads} --topPercent 0.5 \
                        --lcaAlgorithm longReads --longReads true --minSupportPercent 0 --verbose" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running daa-meganizer" >> ${job}
        echo -e "               touch ${statuspath}/stage17b.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage17b.finished" >> ${job}
        echo -e "       fi" >> ${job}
       	echo -e "	paste <(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_dragonflye_hybrid_contigs_blastx.daa -c2c Taxonomy | awk '{print \$1}') \
				<(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_dragonflye_hybrid_contigs_blastx.daa -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') \
				> ${path}/blast/${sample}_dragonflye_hybrid_contigs_blastx_daa_summary_count.tsv" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage17b=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage17b}")
fi

# STAGE 17c: meganize blastx dragonflye contig output (dependency STAGE 15)
if containsElement "17c" "${stages[@]}";
then
        echo "STAGE 17c: meganize blastx dragonflye LR contig output"
        [ -d ${path}/blast ] && echo "${path}/blast \
                exists. ! creating." || mkdir -p ${path}/blast
        job=${jobpath}/stage17c.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage17c" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --chdir=${megandb}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage17c.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage17c.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage12c} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage12c}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage17c.finished ]]; then" >> ${job}
        echo -e "       ${meganpath}/daa-meganizer --in ${path}/blast/${sample}_dragonflye_LR_contigs_blastx.daa \
			--mapDB ${megandb}/megan-map.db --acc2taxa ${megandb}/ncbi.map --threads ${threads} --topPercent 0.5 \
                        --lcaAlgorithm longReads --longReads true --minSupportPercent 0 --verbose" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running daa-meganizer" >> ${job}
        echo -e "               touch ${statuspath}/stage17c.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage17c.finished" >> ${job}
        echo -e "       fi" >> ${job}
       	echo -e "	paste <(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_dragonflye_LR_contigs_blastx.daa -c2c Taxonomy | awk '{print \$1}') \
				<(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_dragonflye_LR_contigs_blastx.daa -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') \
				> ${path}/blast/${sample}_dragonflye_LR_contigs_blastx_daa_summary_count.tsv" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage17c=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage17c}")
fi

# STAGE 17d: meganize blastx unicycler contig output (dependency STAGE 15)
if containsElement "17d" "${stages[@]}";
then
        echo "STAGE 17d: meganize blastx unicycler LR contig output"
        [ -d ${path}/blast ] && echo "${path}/blast \
                exists. ! creating." || mkdir -p ${path}/blast
        job=${jobpath}/stage17d.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage17d" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --chdir=${megandb}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage17d.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage17d.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage12d} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage12d}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage17d.finished ]]; then" >> ${job}
        echo -e "       ${meganpath}/daa-meganizer --in ${path}/blast/${sample}_unicycler_contigs_blastx.daa \
			--mapDB ${megandb}/megan-map.db --acc2taxa ${megandb}/ncbi.map --threads ${threads} --topPercent 0.5 \
                        --lcaAlgorithm longReads --longReads true --minSupportPercent 0 --verbose" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running daa-meganizer" >> ${job}
        echo -e "               touch ${statuspath}/stage17d.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage17d.finished" >> ${job}
        echo -e "       fi" >> ${job}
       	echo -e "	paste <(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_unicycler_contigs_blastx.daa -c2c Taxonomy | awk '{print \$1}') \
				<(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_unicycler_contigs_blastx.daa -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') \
				> ${path}/blast/${sample}_unicycler_contigs_blastx_daa_summary_count.tsv" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage17d=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage17d}")
fi

# STAGE 18: meganize blastx SPAdes contig output (dependency STAGE 13)
if containsElement "18" "${stages[@]}";
then
        echo "STAGE 18: meganize blastx SPAdes contig output"
        [ -d ${path}/blast ] && echo "${path}/blast \
                exists. ! creating." || mkdir -p ${path}/blast
        job=${jobpath}/stage18.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage18" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --chdir=${megandb}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage18.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage18.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage13} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage13}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage18.finished ]]; then" >> ${job}
        echo -e "       ${meganpath}/daa-meganizer --in ${path}/blast/${sample}_spades_contigs_blastx.daa \
			--mapDB ${megandb}/megan-map.db --acc2taxa ${megandb}/ncbi.map --threads ${threads} --topPercent 0.5 \
                        --minSupportPercent 0 --lcaAlgorithm longReads --longReads true --verbose" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running daa-meganizer" >> ${job}
        echo -e "               touch ${statuspath}/stage18.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage18.finished" >> ${job}
        echo -e "       fi" >> ${job}
       	echo -e "	paste <(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_spades_contigs_blastx.daa -c2c Taxonomy | awk '{print \$1}') \
				<(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_spades_contigs_blastx.daa -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') \
				> ${path}/blast/${sample}_spades_contigs_blastx_daa_summary_count.tsv" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage18=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage18}")
fi

# STAGE 19: meganize blastx read output (dependency STAGE 14)
if containsElement "19" "${stages[@]}";
then
       	echo "STAGE 19: meganize blastx read output"
       	[ -d ${path}/blast ] && echo "${path}/blast \
		exists. ! creating." || mkdir -p ${path}/blast
	job=${jobpath}/stage19.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage19" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --chdir=${megandb}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage19.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage19.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage14} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage14}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage19.finished ]]; then" >> ${job}
       	echo -e "	${meganpath}/daa-meganizer --in ${path}/blast/${sample}_reads_blastx.daa --only Taxonomy \
			--mapDB ${megandb}/megan-map.db --acc2taxa ${megandb}/ncbi.map --threads ${threads} --minSupportPercent 0 --topPercent 0.5 \
			--lcaAlgorithm weighted --longReads false --verbose" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running daa-meganizer" >> ${job}
       	echo -e "		touch ${statuspath}/stage19.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage19.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	paste <(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_reads_blastx.daa -c2c Taxonomy | awk '{print \$1}') \
				<(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_reads_blastx.daa -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') \
				> ${path}/blast/${sample}_reads_blastx_daa_summary_count.tsv" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage19=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage19}")
fi

# STAGE 19b: meganize blastx long read output (dependency STAGE 14)
if containsElement "19b" "${stages[@]}";
then
       	echo "STAGE 19b: meganize blastx read output"
       	[ -d ${path}/blast ] && echo "${path}/blast \
		exists. ! creating." || mkdir -p ${path}/blast
	job=${jobpath}/stage19b.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage19b" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --chdir=${megandb}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage19b.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage19b.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage14b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage14b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage19b.finished ]]; then" >> ${job}
       	echo -e "	${meganpath}/daa-meganizer --in ${path}/blast/${sample}_LR_blastx.daa --only Taxonomy \
			--mapDB ${megandb}/megan-map.db --acc2taxa ${megandb}/ncbi.map --threads ${threads} --minSupportPercent 0 --topPercent 0.5 \
			--lcaAlgorithm weighted --longReads false --verbose" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running daa-meganizer" >> ${job}
       	echo -e "		touch ${statuspath}/stage19b.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage19b.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	paste <(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_LR_blastx.daa -c2c Taxonomy | awk '{print \$1}') \
				<(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_LR_blastx.daa -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') \
				> ${path}/blast/${sample}_LR_blastx_daa_summary_count.tsv" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage19b=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage19b}")
fi

# STAGE 20: meganize blastx metaviral SPAdes contig output (dependency STAGE 15)
if containsElement "20" "${stages[@]}";
then
        echo "STAGE 20: meganize blastx metaviral SPAdes contig output"
        [ -d ${path}/blast ] && echo "${path}/blast \
                exists. ! creating." || mkdir -p ${path}/blast
        job=${jobpath}/stage20.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage20" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}g" >> ${job}
	echo -e "#SBATCH --chdir=${megandb}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage20.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage20.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage15} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage15}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage20.finished ]]; then" >> ${job}
        echo -e "       ${meganpath}/daa-meganizer --in ${path}/blast/${sample}_metaviral_contigs_blastx.daa \
			--mapDB ${megandb}/megan-map.db --acc2taxa ${megandb}/ncbi.map --threads ${threads} --topPercent 0.5 \
                        --lcaAlgorithm longReads --longReads true --minSupportPercent 0 --verbose" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running daa-meganizer" >> ${job}
        echo -e "               touch ${statuspath}/stage20.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage20.finished" >> ${job}
        echo -e "       fi" >> ${job}
       	echo -e "	paste <(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_metaviral_contigs_blastx.daa -c2c Taxonomy | awk '{print \$1}') \
				<(${meganpath}/daa2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_metaviral_contigs_blastx.daa -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\t' OFS='\t') \
				> ${path}/blast/${sample}_metaviral_contigs_blastx_daa_summary_count.tsv" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage20=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage20}")
fi

# STAGE 21: meganize blastn contig output (dependency STAGE 16)
if containsElement "21" "${stages[@]}";
then
       	echo "STAGE 21: meganize blastn contig output"
       	[ -d ${path}/blast ] && echo "${path}/blast \
		exists. ! creating." || mkdir -p ${path}/blast
	job=${jobpath}/stage21.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage21" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --mem=${memory}" >> ${job}
	echo -e "#SBATCH --chdir=${megandb}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage21.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage21.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage16} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage16}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage21.finished ]]; then" >> ${job}
       	echo -e "	${meganpath}/blast2rma --in ${path}/blast/${sample}_metaspades_contigs.mmseqs.out \
			--reads ${path}/spades/meta_pe_trim/contigs.fasta \
			--acc2taxa ${megandb}/ncbi.map --topPercent 0.5 \
			--longReads true --format BlastTab --blastMode BlastN --lcaAlgorithm longReads \
			--mapDB ${megandb}/megan-nucl.db --threads ${threads} --minSupportPercent 0 \
			--out ${path}/blast/${sample}_metaspades_contigs_blastn.rma6 --verbose" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running blast2rma" >> ${job}
       	echo -e "		touch ${statuspath}/stage21.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage21.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "	paste <(${meganpath}/rma2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_metaspades_contigs_blastn.rma6 -c2c Taxonomy | awk '{print \$1}') \
				<(${meganpath}/rma2info -P ${meganpath}/.MEGAN.def -i ${path}/blast/${sample}_metaspades_contigs_blastn.rma6 -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\t' OFS='\t') \
				> ${path}/blast/${sample}_metaspades_contigs_blastn_rma6_summary_count.tsv" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage21=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage21}")
fi

# STAGE 22: metaQUAST metaSPAdes assembly validation (dependency STAGE 6)
if containsElement "22" "${stages[@]}";
then
       	echo "STAGE 22: metaQUAST metaSPAdes assembly validation"
       	[ -d ${path}/quast/metaquast ] && \
		echo "${path}/quast/metaquast \
		exists. ! creating." || mkdir -p \
		${path}/quast/metaquast
	job=${jobpath}/stage22.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage22" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage22.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage22.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage6} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage6}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage22.finished ]]; then" >> ${job}
       	echo -e "	metaquast.py -t ${threads} --circos \
			--rna-finding --circos --mgm --conserved-genes-finding \
			--pe12 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
			-o ${path}/quast/metaquast \
		       	${path}/spades/meta_pe_trim/contigs.fasta" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running metaquast" >> ${job}
       	echo -e "		touch ${statuspath}/stage22.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage22.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage22=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage22}")
fi

# STAGE 23: QUAST SPAdes assembly validation (dependency STAGE 7)
if containsElement "23" "${stages[@]}";
then
        echo "STAGE 23: QUAST SPAdes assembly validation"
        [ -d ${path}/quast/quast ] && \
                echo "${path}/quast/quast \
                exists. ! creating." || mkdir -p \
                ${path}/quast/quast
        job=${jobpath}/stage23.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage23" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage23.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage23.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage7} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage7}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage23.finished ]]; then" >> ${job}
        echo -e "       quast.py -t ${threads} --circos \
                        --rna-finding --circos --glimmer --conserved-genes-finding \
                        --pe12 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
                        -o ${path}/quast/quast \
                        ${path}/spades/pe_trim/contigs.fasta" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running quast" >> ${job}
        echo -e "               touch ${statuspath}/stage23.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage23.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage23=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage23}")
fi

# STAGE 23: QUAST dragonflye hybrid assembly validation (dependency STAGE 7)
if containsElement "23b" "${stages[@]}";
then
        echo "STAGE 23b: metaQUAST dragonflye hybrid assembly validation"
        [ -d ${path}/quast/dragonflye_hybrid_metaquast ] && \
                echo "${path}/quast/dragonflye_hybrid_metaquast \
                exists. ! creating." || mkdir -p \
                ${path}/quast/dragonflye_hybrid_metaquast
        job=${jobpath}/stage23b.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage23b" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage23b.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage23b.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage6b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage6b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage23b.finished ]]; then" >> ${job}
        echo -e "       metaquast.py -t ${threads} --circos \
                        --rna-finding --circos --glimmer --conserved-genes-finding \
                        --pe12 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
                        --nanopore ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.fastq \
                        -o ${path}/quast/dragonflye_hybrid_metaquast \
                        ${path}/dragonflye/hybrid/contigs.fa" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running quast" >> ${job}
        echo -e "               touch ${statuspath}/stage23b.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage23b.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage23b=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage23b}")
fi

# STAGE 23: QUAST dragonflye hybrid assembly validation (dependency STAGE 7)
if containsElement "23bb" "${stages[@]}";
then
        echo "STAGE 23bb: QUAST dragonflye hybrid assembly validation"
        [ -d ${path}/quast/dragonflye_hybrid_quast ] && \
                echo "${path}/quast/dragonflye_hybrid_quast \
                exists. ! creating." || mkdir -p \
                ${path}/quast/dragonflye_hybrid_quast
        job=${jobpath}/stage23bb.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage23bb" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage23bb.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage23bb.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage6b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage6b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage23bb.finished ]]; then" >> ${job}
        echo -e "       quast.py -t ${threads} --circos \
                        --rna-finding --circos --glimmer --conserved-genes-finding \
                        --pe12 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
                        --nanopore ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.fastq \
                        -o ${path}/quast/dragonflye_hybrid_quast \
                        ${path}/dragonflye/hybrid/contigs.fa" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running quast" >> ${job}
        echo -e "               touch ${statuspath}/stage23bb.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage23bb.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage23bb=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage23bb}")
fi

# STAGE 23: QUAST dragonflye LR assembly validation (dependency STAGE 7)
if containsElement "23c" "${stages[@]}";
then
        echo "STAGE 23c: metaQUAST dragonflye LR assembly validation"
        [ -d ${path}/quast/dragonflye_LR_metaquast ] && \
                echo "${path}/quast/dragonflye_LR_metaquast \
                exists. ! creating." || mkdir -p \
                ${path}/quast/dragonflye_LR_metaquast
        job=${jobpath}/stage23c.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage23c" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage23c.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage23c.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage6b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage6b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage23c.finished ]]; then" >> ${job}
        echo -e "       metaquast.py -t ${threads} --circos \
                        --rna-finding --circos --glimmer --conserved-genes-finding \
                        --nanopore ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.fastq \
                        -o ${path}/quast/dragonflye_LR_metaquast \
                        ${path}/dragonflye/LR/contigs.fa" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running quast" >> ${job}
        echo -e "               touch ${statuspath}/stage23c.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage23c.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage23c=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage23c}")
fi

# STAGE 23: QUAST dragonflye LR assembly validation (dependency STAGE 7)
if containsElement "23cc" "${stages[@]}";
then
        echo "STAGE 23cc: QUAST dragonflye LR assembly validation"
        [ -d ${path}/quast/dragonflye_LR_quast ] && \
                echo "${path}/quast/dragonflye_LR_quast \
                exists. ! creating." || mkdir -p \
                ${path}/quast/dragonflye_LR_quast
        job=${jobpath}/stage23cc.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage23cc" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage23cc.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage23cc.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage6b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage6b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage23cc.finished ]]; then" >> ${job}
        echo -e "       quast.py -t ${threads} --circos \
                        --rna-finding --circos --glimmer --conserved-genes-finding \
                        --nanopore ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.fastq \
                        -o ${path}/quast/dragonflye_LR_quast \
                        ${path}/dragonflye/LR/contigs.fasta" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running quast" >> ${job}
        echo -e "               touch ${statuspath}/stage23cc.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage23cc.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage23cc=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage23cc}")
fi

# STAGE 23d: QUAST unicycler assembly validation (dependency STAGE 7)
if containsElement "23d" "${stages[@]}";
then
        echo "STAGE 23: QUAST unicycler hybrid isolate assembly validation"
        [ -d ${path}/quast/unicycler_hybrid_quast ] && \
                echo "${path}/quast/unicycler_hybrid_quast \
                exists. ! creating." || mkdir -p \
                ${path}/quast/unicycler_hybrid_quast
        job=${jobpath}/stage23d.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage23d" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage23d.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage23d.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage6b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage6b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage23d.finished ]]; then" >> ${job}
        echo -e "       quast.py -t ${threads} --circos \
                        --rna-finding --circos --glimmer --conserved-genes-finding \
                        --pe12 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
                        --nanopore ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_LR.fastq \
                        -o ${path}/quast/unicycler_hybrid_quast \
                        ${path}/unicycler/isolate/assembly.fasta" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running quast" >> ${job}
        echo -e "               touch ${statuspath}/stage23d.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage23d.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage23d=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage23d}")
fi

# STAGE 24: QUAST metaviral SPAdes assembly validation (dependency STAGE 8)
if containsElement "24" "${stages[@]}";
then
        echo "STAGE 24: QUAST SPAdes assembly validation"
        [ -d ${path}/quast/metaviral ] && \
                echo "${path}/quast/metaviral \
                exists. ! creating." || mkdir -p \
                ${path}/quast/metaviral
        job=${jobpath}/stage23.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage24" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage24.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage24.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage8} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage8}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
        echo -e "if [[ ! -e ${statuspath}/stage24.finished ]]; then" >> ${job}
        echo -e "       quast.py -t ${threads} --circos \
                        --rna-finding --circos --gene-finding --conserved-genes-finding \
                        --pe12 ${path}/trim/shotgun/${sample}_host_contaminant_rRNA_removed_pe.fastq.gz \
                        -o ${path}/quast/metaviral \
                        ${path}/spades/metaviral_pe_trim/contigs.fasta" >> ${job}
        echo -e "       out=\$?" >> ${job}
        echo -e "       if [[ \${out} -ne 0 ]]; then" >> ${job}
        echo -e "               echo error \${out} running quast" >> ${job}
        echo -e "               touch ${statuspath}/stage24.not_finished" >> ${job}
        echo -e "       else" >> ${job}
        echo -e "               touch ${statuspath}/stage24.finished" >> ${job}
        echo -e "       fi" >> ${job}
        echo -e "fi" >> ${job}
        echo "submitting ${job} to SLURM queue"
	jobid_stage24=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage24}")
fi

# STAGE 25: produce posttrim fastqc reports (dependency STAGE 5)
if containsElement "25" "${stages[@]}";
then
       	echo "STAGE 25: producing FastQC posttrim reports"
       	[ -d ${path}/fastqc/posttrim ] && \
		echo "${path}/fastqc/posttrim/shotgun exists. ! creating." \
		|| mkdir -p ${path}/fastqc/posttrim
	job=${jobpath}/stage25.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage25" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage25.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage25.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	if [[ ${jobid_stage5} && ${jobid_stage5b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage5}:${jobid_stage5b}" >> ${job}
	elif [[ ${jobid_stage5} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage5}" >> ${job}
	elif [[ ${jobid_stage5b} ]]; then
	       	echo -e "#SBATCH --dependency=afterok:${jobid_stage5b}" >> ${job}
	fi

	echo -e "eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"" >> ${job}
        echo -e "conda activate md" >> ${job}

        echo -e "set -e" >> ${job}
	echo -e "if [[ ! -e ${statuspath}/stage25.finished ]]; then" >> ${job}
       	echo -e "	fastqc --outdir ${path}/fastqc/posttrim \
			--java java \
			--threads ${threads} ${path}/trim/shotgun/*.fastq*" >> ${job}
       	echo -e "	out=\$?" >> ${job}
       	echo -e	"	if [[ \${out} -ne 0 ]]; then" >> ${job}
       	echo -e "		echo error \${out} running fastqc" >> ${job}
       	echo -e "		touch ${statuspath}/stage25.not_finished" >> ${job}
       	echo -e "	else" >> ${job}
       	echo -e "		touch ${statuspath}/stage25.finished" >> ${job}
       	echo -e "	fi" >> ${job}
       	echo -e "fi" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage25=$(sbatch ${job} | awk '{ print $NF }')
	jobscript+=("${jobid_stage25}")
fi

# STAGE 26: test job completion
if containsElement "26" "${stages[@]}";
then
       	echo "STAGE 26: test job completion"
	printf -v jobscriptcsv '%s,' "${jobscript[@]}"
	job=${jobpath}/stage26.sh
	echo -e "#!/bin/bash" > ${job}
	echo -e "#SBATCH --job-name=MD:${sample_path}:stage26" >> ${job}
	echo -e "#SBATCH --partition=${partition}" >> ${job}
	echo -e "#SBATCH --account=${account}" >> ${job}
	echo -e "#SBATCH --cpus-per-task=${threads}" >> ${job}
	echo -e "#SBATCH --output=${logpath}/stage26.out" >> ${job}
	echo -e "#SBATCH --error=${logpath}/stage26.err" >> ${job}
	echo -e "#SBATCH --open-mode=append" >> ${job}
	echo -e "#SBATCH --dependency=afterany:${jobscriptcsv%,}" >> ${job}
        echo -e "set -e" >> ${job}
	for num in ${stages[@]}
	do
	       	[[ ${num} != ${last_stage} ]] && echo -e "[[ -e ${statuspath}/stage${num}.finished ]] && \\" >> ${job}
	done
	echo -e "echo ${sample} assembly pipeline completed > ${statuspath}/assembly_pipeline.finished || \\" >> ${job}
	echo -e "echo ${sample} did not complete assembly pipeline > ${statuspath}/assembly_pipeline.finished" >> ${job}
	echo "submitting ${job} to SLURM queue"
	jobid_stage26=$(sbatch ${job} | awk '{ print $NF }')
fi

echo "Finishing process for ${sample} at $today"
echo "*****************************"
exit 1
