#!/bin/bash
#bash /home/scripts/GenerateTxpFastaFromDecoy.sh N2_mt-Co2-201.txpnames N2_mt-Co2-201

decoy_inputfile=$1
subsampletxp_outputfile=$2
while IFS= read -r line
do
	sed -n "/>${line}/,/^>/p" /home/0309meeting/assist_issac/from_issac/genomeDir/gencode.vM25.transcripts.fa >> ${subsampletxp_outputfile}.fa
	sed -i "$ d" ${subsampletxp_outputfile}.fa
done < "$decoy_inputfile"
