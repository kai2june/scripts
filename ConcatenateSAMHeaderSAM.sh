#!/bin/bash
#bash /home/scripts/ConcatenateSAMHeaderSAM.sh N2_mt-Co2-201.samheader N2_mt-Co2-201.noheadersam N2_mt-Co2-201

samheader_inputfile=$1
sam_inputfile=$2
sam_outputfile=$3

cat ${samheader_inputfile} ${sam_inputfile} > ${sam_outputfile}.sam
