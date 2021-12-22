#!/bin/bash
# author : jimmy zhuang
# date 2021/07/04

percent=0
basepath="/home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep"
tool="jimmy_salmon"
index="txp.add_nascent.multimapping.seqBias.gcBias.1221"
# branch="transcript_in_N_eqvclass"
mode="idx"

replicate=10
jimmy=""
for x in $(seq 1 ${replicate})
do
    echo "${basepath}${x}/${tool}.${index}/${tool}.${index}.${mode}quant/logs/salmon_quant.log"
	jimmy=${jimmy}+`grep 'Mapping rate' ${basepath}${x}/${tool}.${index}/${tool}.${index}.${mode}quant/logs/salmon_quant.log | awk '{printf($8" ")}'`
done

echo ${jimmy} | awk -F"% " '{for(i=1; i<NF+1; ++i) {print $i;sum+=$i;}}END{printf ("avg:"sum/NF"\n")}'

