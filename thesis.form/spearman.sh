#!/bin/bash

percent=0
basepath="/home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep"
tool="jimmy_salmon"
index="txp.add_nascent.multimapping.seqBias.gcBias.1221"
TPMdir="computeTPM" # "computeTPM.txp" means only txp for correlation, it's for nascent0% non-baseline experiments; "computeTPM" means txp + nascent for correlation
# tool="jimmy_salmon"
# index="txp.add_nascent.all_improvements.seqBias.gcBias"
# branch="transcript_in_N_eqvclass"
mode="idx"

replicate=10
jimmy=""
for x in $(seq 1 ${replicate})
do
    if [ ${tool} = "rsem" ]
        then
            echo "${basepath}${x}/${tool}.${index}/${tool}.${index}.quant.stat/${TPMdir}/correlation.txt"
            jimmy=${jimmy}+`awk 'NR==1' ${basepath}${x}/${tool}.${index}/${tool}.${index}.quant.stat/${TPMdir}/correlation.txt | awk '{printf($2" ")}'`
    elif [ ${tool} = "express" ] || [ ${tool} = "cufflinks" ] || [ ${tool} = "sailfish" ] || [ ${tool} = "kallisto" ]
        then
            echo "${basepath}${x}/${tool}.${index}/${tool}.${index}.quant/${TPMdir}/correlation.txt"
            jimmy=${jimmy}+`awk 'NR==1' ${basepath}${x}/${tool}.${index}/${tool}.${index}.quant/${TPMdir}/correlation.txt | awk '{printf($2" ")}'`
    elif [ ${tool} = "jimmy_salmon" -o ${tool} = "salmon" ]
        then
            echo "${basepath}${x}/${tool}.${index}/${tool}.${index}.${mode}quant/${TPMdir}/correlation.txt"
            jimmy=${jimmy}+`awk 'NR==1' ${basepath}${x}/${tool}.${index}/${tool}.${index}.${mode}quant/${TPMdir}/correlation.txt | awk '{printf($2" ")}'`
    fi
done

echo ${jimmy} | awk -F"+" '{for(i=1; i<NF+1; ++i) {print $i;sum+=$i;}}END{printf("avg:%.4f\n",sum/(NF-1))}'