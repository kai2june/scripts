#!/bin/bash
# author : jimmy zhuang
# date 2021/07/04

basepath="/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep"
tool="jimmy_salmon"
index="txp"
branch="master"
mode="idx"

replicate=10
jimmy=""
for x in $(seq 1 ${replicate})
do
	jimmy=${jimmy}+`grep 'Mapping rate' ${basepath}${x}/${tool}.${index}/${tool}.${index}.${branch}.${mode}quant/logs/salmon_quant.log | awk '{printf($8" ")}'`
done

echo ${jimmy} | awk -F"% " '{for(i=1; i<NF+1; ++i) {print $i;sum+=$i;}}END{printf ("avg:"sum/NF"\n")}'

