#!/bin/bash
#bash /home/scripts/GenerateSAMHeaderFromDecoy.sh N2_mt-Co2-201.txpnames ../N1.sam N2_mt-Co2-201

decoy_inputfile=$1
sam_inputfile=$2
samheader_outputfile=$3

head -n 1 ${sam_inputfile} | grep "^@HD" >> ${samheader_outputfile}.samheader
while IFS= read -r line
do
	grep -m 1 ${line} ${sam_inputfile} >> ${samheader_outputfile}.samheader
done < ${decoy_inputfile}
