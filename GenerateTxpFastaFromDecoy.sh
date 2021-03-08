#!/bin/bash
###bash /home/scripts/GenerateTxpFastaFromDecoy.sh /home/0309meeting/assist_issac/from_issac/hisat2/transcriptome_N2_sortreadname_output.txpnames mytmp.fa

decoy_inputfile=$1
subsampletxp_outputfile=$2
while IFS= read -r line
do
	sed -n "/>${line}/,/^>/p" /home/0309meeting/assist_issac/from_issac/genomeDir/gencode.vM25.transcripts.fa >> ${subsampletxp_outputfile}
	sed -i "$ d" ${subsampletxp_outputfile}
done < "$decoy_inputfile"
