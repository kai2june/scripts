#!/bin/bash

replicate=10
percent=0
tool="express"
index="txp.salmonsam.EM100"
for i in $(seq 1 ${replicate})
do
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i} && 
    if [ ! -d ${tool}.${index} ] 
        then
            mkdir ${tool}.${index}
    fi
    cd ${tool}.${index} &&
    express --fr-stranded -O 100 --output-dir ${tool}.${index}.quant /home/droso.reference/gffread.txp.fa ../salmon.txp/salmon.txp.idxquant.sam &&
    cd ${tool}.${index}.quant &&
    rm -rf computeTPM &&
    mkdir computeTPM &&
    cd computeTPM &&
    python3 /home/scripts/ComputeTPMCorrelation.py -g /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%.pro -b /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%.bed -e ../results.xprs & 
done



