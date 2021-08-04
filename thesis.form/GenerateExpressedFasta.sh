#!/bin/bash

for i in $(seq 1 10)
do
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep${i} &&
    # awk '{if($10 > 0) print $2}' txpgene0%.pro > expressed.txt &&
    # sed -i 's/^/transcript_id "/; s/$/"/' expressed.txt &&
    # grep -f expressed.txt /home/droso.reference/Drosophila_melanogaster.BDGP6.80_txp.gtf > expressed.gtf &&
    # gffread -w expressed.fa -g /home/droso.reference/Drosophila_melanogaster.BDGP6.dna.chromosome.all.fa expressed.gtf &&
    # salmon index -t expressed.fa -i salmon.txp.expressed.index -p 60 --keepDuplicates &&
    # mkdir salmon.txp.expressed && mv expressed.fa expressed.txt expressed.gtf salmon.txp.expressed.index salmon.txp.expressed/ && cd salmon.txp.expressed && 
    # salmon quant -i salmon.txp.expressed.index -l A -1 ../txpgene0%_1.fasta -2 ../txpgene0%_2.fasta -d -o salmon.txp.expressed.idxquant &&
    cd salmon.txp.expressed/salmon.txp.expressed.idxquant && rm -rf computeTPM && mkdir computeTPM && cd computeTPM &&
    python3 /home/scripts/ComputeTPMCorrelation.py -g ../../../txpgene0%.pro -b ../../../txpgene0%.bed -t ../quant.sf &
done