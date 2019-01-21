#!/bin/bash
sample_path=$1
out_path=$2
output=${out_path}/QC_calculate.txt

printf "total_reads,align_success,align_overlap,TRA,TRB\n" >> ${output}

for sample_dir in `ls -l ${sample_path} | awk '/^d/{print $NF}'`
do
	align_report=${sample_path}/${sample_dir}/alignmentReport.txt
	total_reads=`cat ${align_report} | grep "Total sequencing reads" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	align_success=`cat ${align_report} | grep "Successfully aligned reads" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	align_overlap=`cat ${align_report} | grep "Overlapped and aligned" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	TRA=`cat ${align_report} | grep "TRA chains" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	TRB=`cat ${align_report} | grep "TRB chains" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	printf "%s,%s,%s,%s,%s\n" "${total_reads}" "${align_success}" "${align_overlap}" "${TRA}" "${TRB}" >> ${output}
done
