#!/bin/bash
source md.config
basedir=${md_datadir}
shopt -s extglob

timestamp() { 
	date +%Y-%m-%d_%H%M%S 
}

display_usage() { 
	echo -e "\nUsage:\n$0 [-a assembly_flag] [-b bbmap options] [-d host_db for host removal] [-h help] [-j assembly_type] [-k bbduk options] [-m memory in GB] [-o output_directory] [-r readlist] [-s stages] [-t threads]\n" >&2
}

while getopts "h?a:b:d:k:m:o:r:s:t:j:" opt; do
       	case "$opt" in
		a)	assembly_flag="$OPTARG" # ("a" for all assembly stages,"i" for isolate assembly,"m" for metagenomic assembly,"t" for isolate with targeted reference matches,"v" for metaviral assembly)
			;;
		b)	bbmap_opts="${OPTARG@Q}" # (optional arguments for bbmap)
			;;
		d)	host_db="$OPTARG" # (path/to/host/fasta file for removal or targeted reference matching)
			;;
	       	h|\?)	display_usage 
			exit 1
			;;
		k)	bbduk_opts="${OPTARG@Q}" # (optional arguments for bbduk)
			;;
	       	m)  	memory="$OPTARG" ;; # memory per STAGE in GB e.g. 128
	       	o)  	outdir="$basedir/MD_output/$OPTARG" ;; # e.g. batch1
	       	r)  	readlist="$basedir/readlist/$OPTARG" ;; # file which contains path-to-reads
	       	s)  	csv_stages="$OPTARG" ;; # comma separated list of stages to run; null for all
	       	t)  	threads="$OPTARG" ;; # threads per stage e.g. 64
	       	j)	assembly_type="$OPTARG" ;; # type of assembly according to reads e.g. s (short), h (hybrid), l (long)
       	esac
done

if ((OPTIND == 1))
then 
	echo "No options specified" >&2
	display_usage
	exit 1
fi

shift $((OPTIND-1))

if (($# == 1))
then
	echo "No positional arguments required" >&2
	display_usage
	exit 1
fi

if [[ -z ${readlist} ]]
then
	echo "-r readlist required" >&2
	exit 1
fi

[[ ! -z ${host_db} ]] && host_db_flag="-d ${host_db}"
[[ ! -z ${bbduk_opts} ]] && bbduk_opts_flag="-k ${bbduk_opts}"
[[ ! -z ${bbmap_opts} ]] && bbmap_opts_flag="-b ${bbmap_opts}"
[[ -z ${outdir} ]] && outdir=$basedir/MD_output/out.$(timestamp)
[[ ! -z ${assembly_flag} ]] && aflag="-a ${assembly_flag}" || aflag="-a m"
[[ ! -z ${assembly_type} ]] && atype="-j ${assembly_type}" || atype="-j s"
[[ ! -z ${csv_stages} ]] && stages="-s ${csv_stages}"
[[ -z ${memory} ]] && memory=128
[[ -z ${threads} ]] && threads=62
[[ -z ${assembly_type} ]] && assembly_type="s"

host_db_flag=$(echo ${host_db_flag} | xargs) # remove whitespace
assembly_flag=$(echo ${assembly_flag} | xargs) # remove whitespace
assembly_type=$(echo ${assembly_type} | xargs) # remove whitespace
bbmap_opts_flag=$(echo ${bbmap_opts_flag} | xargs) # remove whitespace
stages=$(echo ${stages} | xargs) # remove whitespace
mkdir -p ${outdir}
log="$outdir/queue.log"
exec >>${log} 2>&1

echo -e "\n*** Starting MetaDetector queue process at $(timestamp) ***"

count_samples=()

for sample in $(cat ${readlist})
do
	sample_name=$(basename ${sample} | \
		sed -n 's:_S[0-9]*_L00[0-9]_R[1-2]_001.fastq.gz::p;s:_S[0-9]*_R[1-2]_001.fastq.gz::p;s:_[Rr][1-2].fastq.gz::p;s:_[Rr][1-2]_001.fastq.gz::p;s:_[1-2].fastq.gz::p;s:.fastq.gz::p' )
	[ ! -d ${outdir}/${sample_name} ] && mkdir -p ${outdir}/${sample_name}
	count_samples+="$sample_name "
	echo $count_samples
	if [[ "$assembly_type" == "s" ]]
	then
		r1=(${outdir}/${sample_name}/${sample_name}_*[rR]1*)
		r2=(${outdir}/${sample_name}/${sample_name}_*[rR]2*)
		if [[ ! -s ${r1} ]] || [[ ! -s ${r2} ]]
			then
			echo -e "Copying ${sample} at $(timestamp)\n"
			cp -rvL ${sample} ${outdir}/${sample_name}
		fi
		if [ $(grep -ow $sample_name <<< ${count_samples[@]} | wc -l) = 2 ]
		then
			ln -rsf ${r1} ${outdir}/${sample_name}/${sample_name}_R1.fastq.gz 
			ln -rsf ${r2} ${outdir}/${sample_name}/${sample_name}_R2.fastq.gz 
			echo -e "Queueing ${sample_name} at $(timestamp)\n"
			echo ${basedir}/proc_assembly.sh -r ${outdir}/${sample_name} -m "${memory}" -t "${threads}" "${aflag}" "${host_db_flag}" \
				"${stages}" "${atype}" "${bbduk_opts_flag}" "${bbmap_opts_flag}"
			${basedir}/proc_assembly.sh -r ${outdir}/${sample_name} -m "${memory}" -t "${threads}" "${aflag}" "${atype}" "${host_db_flag}" "${stages}" "${bbduk_opts_flag}" "${bbmap_opts_flag}"
			#echo ${outdir}/${sample_name}
		fi
	elif [[ "$assembly_type" == "l" ]]
	then 
		r=(${outdir}/${sample_name}/${sample_name}.fastq.gz)
		if [[ ! -s ${r} ]]
			then
			echo -e "Copying ${sample} at $(timestamp)\n"
			cp -rvL ${sample} ${outdir}/${sample_name}
		fi

		echo -e "Queueing ${sample_name} at $(timestamp)\n"
		echo ${basedir}/proc_assembly.sh -r ${outdir}/${sample_name} -m "${memory}" -t "${threads}" "${aflag}" "${atype}" "${host_db_flag}" "${bbmap_opts_flag}" "${stages}"
		${basedir}/proc_assembly.sh -r ${outdir}/${sample_name} -m "${memory}" -t "${threads}" "${aflag}" "${atype}" "${host_db_flag}" "${bbmap_opts_flag}" "${stages}"
		#echo ${outdir}/${sample_name}
	fi
	if [[ "$assembly_type" == "h" ]]
	then
		r=(${outdir}/${sample_name}/${sample_name}.fastq.gz)
		r1=(${outdir}/${sample_name}/${sample_name}_*[rR]1*)
		r2=(${outdir}/${sample_name}/${sample_name}_*[rR]2*)
		if [[ ! -s ${r1} ]] || [[ ! -s ${r2} ]] || [[ ! -s ${r} ]]
			then
			echo -e "Copying ${sample} at $(timestamp)\n"
			cp -rvL ${sample} ${outdir}/${sample_name}
		fi
		if [ $(grep -ow $sample_name <<< ${count_samples[@]} | wc -l) = 3 ]
		then
			ln -rsf ${r1} ${outdir}/${sample_name}/${sample_name}_R1.fastq.gz 
			ln -rsf ${r2} ${outdir}/${sample_name}/${sample_name}_R2.fastq.gz 
			
			echo -e "Queueing ${sample_name} at $(timestamp)\n"
			echo ${basedir}/proc_assembly.sh -r ${outdir}/${sample_name} -m "${memory}" -t "${threads}" "${aflag}" "${atype}" "${host_db_flag}" "${bbmap_opts_flag}" "${stages}"
			${basedir}/proc_assembly.sh -r ${outdir}/${sample_name} -m "${memory}" -t "${threads}" "${aflag}" "${atype}" "${host_db_flag}" "${bbmap_opts_flag}" "${stages}"
			#echo ${outdir}/${sample_name}
		fi
	fi	
done
echo -e "Input array = ${count_samples[@]}"
echo -e "*** Finished at $(timestamp) ***"
exit 0
