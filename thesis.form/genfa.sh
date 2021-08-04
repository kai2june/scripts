#!/bin/bash

txp="FBtr0091512"
txp2="FBtr0334671"
grep "transcript_id \"${txp}\"" /home/droso.reference/Drosophila_melanogaster.BDGP6.80_txp.gtf > ${txp}.gtf && 
gffread -w ${txp}.fa -g /home/droso.reference/Drosophila_melanogaster.BDGP6.dna.chromosome.all.fa ${txp}.gtf && 
rm ${txp}.gtf
grep "transcript_id \"${txp2}\"" /home/droso.reference/Drosophila_melanogaster.BDGP6.80_txp.gtf > ${txp2}.gtf && 
gffread -w ${txp2}.fa -g /home/droso.reference/Drosophila_melanogaster.BDGP6.dna.chromosome.all.fa ${txp2}.gtf && 
rm ${txp2}.gtf
vimdiff -R ${txp}.fa ${txp2}.fa