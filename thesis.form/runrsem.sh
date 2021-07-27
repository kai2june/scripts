#!/bin/bash

replicate=10
percent=0
tool="rsem"
index="txp"
mode="aln"
for i in $(seq 1 ${replicate})
do
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i} && 
    if [ ! -d ${tool}.${index}.${mode} ] 
        then
            mkdir ${tool}.${index}.${mode}
    fi
    cd ${tool}.${index}.${mode} && 
    if [ ${mode} = "idx" ]
        then
            # /home/RSEM/rsem-calculate-expression --bowtie2 -p 60 --paired-end /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%_1.fasta /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%_2.fasta --no-qualities /home/droso.reference/${tool}.${index}.index/${tool}.${index}.index ${tool}.${index}.${mode}quant &&
    elif [ ${mode} = "aln" ]
        then
            /home/RSEM/rsem-calculate-expression -p 60 --paired-end --alignments /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/salmon.txp/salmon.txp.idxquant.sam /home/droso.reference/${tool}.${index}.index ${tool}.${index}.${mode}quant
    fi
    cd ${tool}.${index}.quant.stat &&
    rm -rf computeTPM &&
    mkdir computeTPM &&
    cd computeTPM &&
    python3 /home/scripts/ComputeTPMCorrelation.py -g /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%.pro -b /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%.bed -m /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/${tool}.${index}/${tool}.${index}.quant.isoforms.results & 
done
