#!/bin/bash
###bash /home/scripts/ConcatenateSAMHeaderSAM.sh transcriptome_N2_sortreadname_subsample.samheader transcriptome_N2_sortreadname_subsample.sam transcriptome_N2_subsample.sam

samheader_inputfile=$1
sam_inputfile=$2
sam_outputfile=$3

cat ${samheader_inputfile} ${sam_inputfile} > ${sam_outputfile}
