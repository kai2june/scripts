#!/bin/bash

replicate=10
percent=0
tool="kallisto"
index="txp"
for i in $(seq 1 ${replicate})
do
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i} && 
    if [ ! -d ${tool}.${index} ] 
        then
            mkdir ${tool}.${index}
    fi
    cd ${tool}.${index} && 
    kallisto quant -i /home/droso.reference/${tool}.${index}.index --fr-stranded -o ${tool}.${index}.quant /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%_1.fasta /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%_2.fasta &&
    cd ${tool}.${index}.quant &&
    rm -rf computeTPM &&
    mkdir computeTPM &&
    cd computeTPM &&
    python3 /home/scripts/ComputeTPMCorrelation.py -g /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%.pro -b /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%.bed -k ../abundance.tsv & 
done
