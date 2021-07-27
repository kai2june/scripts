#!/bin/bash

transcript=""
while getopts "t:" argv
do
    case "${argv}" in
        t)
            transcript=${OPTARG}
            ;;
        *)
            echo "exception in CreateFasta.sh getopts."
            ;;
    esac
done

grep "transcript_id \"${transcript}\"" ../Drosophila_melanogaster.BDGP6.80_txp.gtf > ${transcript}.gtf
gffread -w ${transcript}.fa -g ../Drosophila_melanogaster.BDGP6.dna.chromosome.all.fa ${transcript}.gtf
rm ${transcript}.gtf
# vimdiff -R FBtr0091882.fa ${transcript}.fa