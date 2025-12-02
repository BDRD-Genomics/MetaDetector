#!/bin/bash
dirc=$1
total=0
cnt=0
for job in $(ls -d ${dirc}/*)
do
	if [[ -f ${job}/log/stage1.out ]] && [[ -f ${job}/status_log/assembly_pipeline.finished ]]; then
		val=$(($(stat -c%Y ${job}/status_log/assembly_pipeline.finished) - $(stat -c%Y ${job}/log/stage1.out)))
		[[ ${val} -lt 0 ]] && val=0
		hours=$(echo "scale=2; ${val}/3600" | bc)
		echo ${job} required ${hours} hours to complete
		total=$(echo ${total} + ${hours} | bc )
		cnt=$((cnt+1))
	fi
done
[[ ${cnt} != 0 ]] && avg=$(echo "scale=2; ${total}/${cnt}" | bc ) && echo ${avg} average hours to completion among ${cnt} jobs
