#!/bin/bash

for i in $(seq 1 10)
do
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene.rep${i}
    mv txp_30M_1.fasta txpgene0%_1.fasta
    mv txp_30M_2.fasta txpgene0%_2.fasta
    mv txp_30M.bed txpgene0%.bed
    mv txp_30M.fasta txpgene0%.fasta
    mv txp_30M.lib txpgene0%.lib
    mv txp_30M.par txpgene0%.par
    mv txp_30M.pro txpgene0%.pro
done