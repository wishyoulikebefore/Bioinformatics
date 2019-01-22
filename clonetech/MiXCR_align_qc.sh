#!/bin/bash

#log_path=/ehpcdata/zuotianyu/clonetech_Jan_15_2019/TCR/SZ_university/align_result
log_path=$1
out_path=$2

output=${out_path}/MiXCR_qc.csv
output_for_analysis=${out_path}/MiXCR_qc_proportion.csv

printf "sample_id,total reads,align success,overlapped and aligned,alignment failed no hits,alignment failed because of absence of V hits,alignment failed because of absence of J hits,no target with both V and J alignments,TRA,TRB\n" >> ${output}
printf "sample_id,align success,overlapped and aligned,alignment failed no hits,alignment failed because of absence of V hits,alignment failed because of absence of J hits,no target with both V and J alignments,TRA,TRB\n" >> ${output_for_analysis}


for sample_id in `ls -l ${log_path} | awk '/^d/{print $NF}'`
do
	align_report=${log_path}/${sample_id}/alignmentReport.txt
	total_reads=`cat ${align_report} | grep "Total sequencing reads" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	align_success=`cat ${align_report} | grep "Successfully aligned reads" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	align_success_rate=`cat ${align_report} | grep "Successfully aligned reads" | cut -f 2 -d ":" | sed 's/[[:space:]]//g' | sed 's/.*(//' | sed 's/%)//'`
	overalp_align=`cat ${align_report} | grep "Overlapped and aligned" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	overalp_align_rate=`cat ${align_report} | grep "Overlapped and aligned" | cut -f 2 -d ":" | sed 's/[[:space:]]//g' | sed 's/.*(//' | sed 's/%)//'`
	align_fail_no_hits=`cat ${align_report} | grep "Alignment failed, no hits" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	align_fail_no_hits_rate=`cat ${align_report} | grep "Alignment failed, no hits" | cut -f 2 -d ":" | sed 's/[[:space:]]//g' | sed 's/.*(//' | sed 's/%)//'`
	align_fail_absence_V=`cat ${align_report} | grep "absence of V hits" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	align_fail_absence_V_rate=`cat ${align_report} | grep "absence of V hits" | cut -f 2 -d ":" | sed 's/[[:space:]]//g' | sed 's/.*(//' | sed 's/%)//'`
	align_fail_absence_J=`cat ${align_report} | grep "absence of J hits" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	align_fail_absence_J_rate=`cat ${align_report} | grep "absence of J hits" | cut -f 2 -d ":" | sed 's/[[:space:]]//g' | sed 's/.*(//' | sed 's/%)//'`
	align_fail_absence_both=`cat ${align_report} | grep "both V and J alignments" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	align_fail_absence_both_rate=`cat ${align_report} | grep "both V and J alignments" | cut -f 2 -d ":" | sed 's/[[:space:]]//g' | sed 's/.*(//' | sed 's/%)//'`
	TRA=`cat ${align_report} | grep "TRA chains" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	TRA_rate=`cat ${align_report} | grep "TRA chains" | cut -f 2 -d ":" | sed 's/[[:space:]]//g' | sed 's/.*(//' | sed 's/%)//'`
	TRB=`cat ${align_report} | grep "TRB chains" | cut -f 2 -d ":" | sed 's/[[:space:]]//g'`
	TRB_rate=`cat ${align_report} | grep "TRB chains" | cut -f 2 -d ":" | sed 's/[[:space:]]//g' | sed 's/.*(//' | sed 's/%)//'`
	printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" "${sample_id}" "${total_reads}" "${align_success}" "${overalp_align}" "${align_fail_no_hits}" "${align_fail_absence_V}" "${align_fail_absence_J}" "${align_fail_absence_both}" "${TRA}" "${TRB}" >> ${output}
	printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" "${sample_id}" "${align_success_rate}" "${overalp_align_rate}" "${align_fail_no_hits_rate}" "${align_fail_absence_V_rate}" "${align_fail_absence_J_rate}" "${align_fail_absence_both_rate}" "${TRA_rate}" "${TRB_rate}" >> ${output_for_analysis}
done
