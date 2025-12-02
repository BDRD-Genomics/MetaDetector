#!/bin/bash
# show stats on assembly pipeline progress

display_usage() {
       	echo -e "\nUsage:\n$0 [-h help] [-d show DependencyNeverSatisfied job ids ] [-f show failed job ids] [-i show incomplete job ids] [-s show completed job ids]" \
	       "<path_to_results>\n" >&2
}

while getopts "dfhsi?" opt; do
       	case "$opt" in
	       	h|\?)
		       	display_usage
		       	exit 0
		       	;;
	       	d)  show_dns=true
		       	;;
	       	f)  show_failed=true
		       	;;
	       	i)  show_inc=true
		       	;;
	       	s)  show_complete=true
		       	;;
       	esac
done

shift $((OPTIND-1))
dir=${1}
account="MD:"$(basename ${dir})

if [ $# -ne 1 ]; then
	display_usage
	exit 0 
fi

levels=$(($(echo $(find ${dir} -type f -path "*/status_log/*" -name "*finished" 2>/dev/null | grep -m 1 -o '/' | uniq -c ) | cut -d' ' -f 1)-1))
done=$(find ${dir} -type f -path "*/status_log/*" -name "assembly_pipeline.finished" -exec \
	grep "assembly pipeline completed" {} + 2>/dev/null | wc -l)
failed_jobs=$(find ${dir} -type f -path "*/status_log/*" -name "assembly_pipeline.finished" -exec \
	grep "did not complete assembly pipeline" {} + 2>/dev/null | wc -l)
dnsjobs=$(squeue -o "%.150j %R" -t PD -A ${account} | tail -n +2 | grep $(basename ${dir}) | grep DependencyNeverSatisfied | cut -d':' -f2 | sort -k1 -u | wc -l)
pending=$(squeue -o "%.150j %R" -t PD -A ${account} | tail -n +2 | grep $(basename ${dir}) | grep -v DependencyNeverSatisfied | cut -d':' -f2 | sort -k1 -u | wc -l)
pending=$(comm -23 <(squeue -o "%.150j %R" -t PD -A ${account} | tail -n +2 | grep $(basename ${dir}) | grep -v DependencyNeverSatisfied | cut -d':' -f2 | sort -k1 -u) \
	<(squeue -o "%.150j %R" -t PD -A ${account} | tail -n +2 | grep $(basename ${dir}) | grep DependencyNeverSatisfied | cut -d':' -f2 | sort -k1 -u) | wc -l)
running=$(squeue -o "%.150j" -t R -A ${account} | tail -n +2 | grep $(basename ${dir}) | cut -d':' -f2 | sort -k1 -u | wc -l)

if [[ ${levels} -gt 0 ]] 
then
     	incjobs=$(comm -3 <(find ${dir} -name "assembly_pipeline.finished" -path "*/status_log/*" -print 2>/dev/null | cut -d'/' -f ${levels} | sort -u) \
	<(find ${dir} -mindepth 1 -maxdepth 1 -type d -not -name "raw_reads" -exec basename {} \; 2>/dev/null| sort -u) | wc -l)
fi

[ -z ${pending} ] && pending=0
[ -z ${running} ] && running=0
[ -z ${dnsjobs} ] && dns=0
[ -z ${incjobs} ] && incjobs=0
[ -z ${failed_jobs} ] && failed_jobs=0

printf "%4s" ${done} " $(basename ${dir}) jobs completed"; echo
printf "%4s" ${failed_jobs} " $(basename ${dir}) jobs finished with failed steps"; echo
printf "%4s" ${dnsjobs} " $(basename ${dir}) jobs DependencyNeverSatisfied"; echo
printf "%4s" ${incjobs} " $(basename ${dir}) jobs incomplete"; echo
printf "%4s" ${running} " $(basename ${dir}) jobs running"; echo
printf "%4s" ${pending} " $(basename ${dir}) jobs pending"; echo

if [ ${show_complete} ]; then
	declare -a compjobs=()
	mapfile -d $'\n' compjobs < <(find ${dir} -type f -name "assembly_pipeline.finished" -exec \
		grep -h "assembly pipeline completed" {} + 2>/dev/null | cut -d' ' -f 1 | sort -h)
	if [ ${#compjobs[@]} -gt 0 ]; then
	       	printf "\n%s\n" "$(basename ${dir}) completed job ids:"
	       	for e in ${compjobs[@]}; do printf "%4s%s\n" '     ' ${e}; done
	fi
fi

if [ ${show_failed} ]; then
	declare -a failed_jobs=()
	mapfile -d $'\n' failed_jobs < \
		<(find ${dir} -type f -name "assembly_pipeline.finished" -exec grep -h "did not complete assembly pipeline" {} + 2>/dev/null | cut -d' ' -f 1 | sort -h)
	if [ ${#failed_jobs[@]} -gt 0 ]; then
	       	printf "\n%s\n" "$(basename ${dir}) failed job ids:"
	       	for e in ${failed_jobs[@]}; do printf "%4s%s\n" '     ' ${e}; done
	fi
fi

if [ ${show_inc} ] && [ ${levels} -gt 0 ]; then
	declare -a incjobs=()
	mapfile -d $'\n' incjobs < <(comm -3 <(find ${dir} -name "assembly_pipeline.finished" -path "*/status_log/*" 2>/dev/null | cut -d'/' -f ${levels} | sort -u) \
	<(find ${dir} -mindepth 1 -maxdepth 1 -type d -not -name "raw_reads" -exec basename {} \; 2>/dev/null | sort -u))
	if [ ${#incjobs[@]} -gt 0 ]; then
	       	printf "\n%s\n" "$(basename ${dir}) incomplete job ids:"
	       	for e in ${incjobs[@]}; do printf "%4s%s\n" '     ' ${e}; done
	fi
fi

if [ ${show_dns} ] && [ ${levels} -gt 0 ]; then
	declare -a dnsjobs=()
	mapfile -d $'\n' dnsjobs < <(squeue -o "%.150j %R" -t PD -A ${account} | tail -n +2 | grep $(basename ${dir}) | grep DependencyNeverSatisfied | sort -k1 -u | awk '{ print $1 }' | cut -d':' -f2 | cut -d'/' -f 2 | sort -u)
	if [ ${#dnsjobs[@]} -gt 0 ]; then
	       	printf "\n%s\n" "$(basename ${dir}) DependencyNeverStatisfied job ids:"
	       	for e in ${dnsjobs[@]}; do printf "%4s%s\n" '     ' ${e}; done
	fi
fi
exit 0
