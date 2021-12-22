#!/bin/bash

replicate=10
percent=60
basepath="/home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep"
index="jimmy_salmon.txp.add_nascent.multimapping.entropy.1220"
TPMdir="computeTPM"
transcript_lower_length=10000
transcript_upper_length=-1

jimmy=""
for x in $(seq 1 ${replicate})
do
    echo "${basepath}${x}/${index}/${index}.idxquant/${TPMdir}/corr_transcript_length_${transcript_lower_length}_${transcript_upper_length}.txt"
    jimmy=${jimmy}+`awk 'NR==2' ${basepath}${x}/${index}/${index}.idxquant/${TPMdir}/corr_transcript_length_${transcript_lower_length}_${transcript_upper_length}.txt | awk '{printf($2" ")}'`
done

echo ${jimmy} | awk -F"+" '{for(i=1; i<NF+1; ++i) {print $i;sum+=$i;}}END{printf("avg:%.4f\n",sum/(NF-1))}'