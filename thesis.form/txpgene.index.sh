#!/bin/bash

for i in $(seq 1 10)
do
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent5%/txpgene5%.rep${i}
    mkdir salmon.txpgene && cd salmon.txpgene
    salmon quant -i /home/droso.reference/salmon.txpgene.index -l A -1 ../txpgene5%_1.fasta -2 ../txpgene5%_2.fasta -d -o salmon.txpgene.quant
done 