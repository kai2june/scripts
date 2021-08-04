#!/bin/bash

for i in $(seq 1 10)
do
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent30% &&
    mkdir txp70%.rep${i} &&
    cd txp70%.rep${i} &&
    cp ../txp70%.par . &&
    flux-simulator -xls -p txp70%.par && 
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent30% &&
    mkdir gene30%.rep${i} &&
    cd gene30%.rep${i} &&
    cp ../gene30%.par . &&
    flux-simulator -xls -p gene30%.par &&
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent30% &&
    mkdir txpgene30%.rep${i} &&
    cd txpgene30%.rep${i} &&
    cat /home/0309meeting/0413/txp_vs_genetxp/nascent30%/txp70%.rep${i}/txp70%.lib /home/0309meeting/0413/txp_vs_genetxp/nascent30%/gene30%.rep${i}/gene30%.lib > txpgene30%.lib &&
    cat /home/0309meeting/0413/txp_vs_genetxp/nascent30%/txp70%.rep${i}/txp70%.pro /home/0309meeting/0413/txp_vs_genetxp/nascent30%/gene30%.rep${i}/gene30%.pro > txpgene30%.pro &&
    cat /home/0309meeting/0413/txp_vs_genetxp/nascent30%/txp70%.rep${i}/txp70%.bed /home/0309meeting/0413/txp_vs_genetxp/nascent30%/gene30%.rep${i}/gene30%.bed > txpgene30%.bed &&
    cat /home/0309meeting/0413/txp_vs_genetxp/nascent30%/txp70%.rep${i}/txp70%.fasta /home/0309meeting/0413/txp_vs_genetxp/nascent30%/gene30%.rep${i}/gene30%.fasta > txpgene30%.fasta &&
    awk 'NR%4==1 || NR%4==2' txpgene30%.fasta > txpgene30%_1.fasta &&
    awk 'NR%4==3 || NR%4==0' txpgene30%.fasta > txpgene30%_2.fasta &
done