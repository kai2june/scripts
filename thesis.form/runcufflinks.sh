#!/bin/bash

replicate=10
percent=0
tool="cufflinks"
index="txp.--dta-cufflinks"
for i in $(seq 2 ${replicate})
do
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i} && 
    if [ ! -f hisat2.genome.--dta-cufflinks.sam ] && [ ! -f hisat2.genome.--dta-cufflinks.sortbySQ.bam ]
        then
            hisat2 -x /home/droso.reference/hisat2.genome.index/hisat2.genome.index -1 txpgene0%_1.fasta -2 txpgene0%_2.fasta -f --dta-cufflinks -p 60 -S hisat2.genome.--dta-cufflinks.sam &&
            samtools sort -t SQ --threads 60 hisat2.genome.--dta-cufflinks.sam > hisat2.genome.--dta-cufflinks.sortbySQ.bam
    fi
    if [ ! -d ${tool}.${index} ] 
        then
            mkdir ${tool}.${index}
    fi
    cd ${tool}.${index} &&
    cufflinks -G /home/droso.reference/Drosophila_melanogaster.BDGP6.80_txp.gtf -p 60 -o ${tool}.${index}.quant ../hisat2.genome.--dta-cufflinks.sortbySQ.bam &&
    cd ${tool}.${index}.quant &&
    rm -rf computeTPM &&
    mkdir computeTPM &&
    cd computeTPM &&
    python3 /home/scripts/ComputeTPMCorrelation.py -g /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%.pro -b /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%.bed -c ../isoforms.fpkm_tracking -i & 
done
