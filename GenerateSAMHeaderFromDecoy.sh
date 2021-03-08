#!/bin/bash
###bash /home/scripts/GenerateSAMHeaderFromDecoy.sh transcriptome_N2_sortreadname_subsample.txpnames transcriptome_N2_sortreadname.sam transcriptome_N2_sortreadname_subsample.samheader

decoy_inputfile=$1
sam_inputfile=$2
samheader_outputfile=$3

head -n 1 ${sam_inputfile} | grep "^@HD" >> ${samheader_outputfile}
while IFS= read -r line
do
	grep -m 1 ${line} ${sam_inputfile} >> ${samheader_outputfile}
done < ${decoy_inputfile}
