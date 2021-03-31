#!/bin/bash
#bash /home/scripts/GenerateTxpFastaFromDecoy.sh N2_mt-Co2-201.txpnames N2_mt-Co2-201

decoy_inputfile=$1
fastainputfile=$2
subsampletxp_outputfile=$3
while IFS= read -r line
do
	sed -n "/>${line}/,/^>/p" ${fastainputfile} >> ${subsampletxp_outputfile}.fa
	sed -i "$ d" ${subsampletxp_outputfile}.fa
done < "$decoy_inputfile"
