#!/bin/bash
#bash /home/scripts/GenerateFastqFromSAM.sh N2_mt-Co2-201.noheadersam N2_reads

sam_inputfile=$1
fastq_outputfile=$2

awk '{print $1, $10, $11}' "$sam_inputfile" | uniq | awk '{printf "@%s\n%s\n+\n%s\n", $1, $2, $3}' > "$fastq_outputfile".fq
awk 'NR%8==1 || NR%8==2 || NR%8==3 || NR%8==4' "$fastq_outputfile".fq > "$fastq_outputfile"_1.fq
awk 'NR%8==5 || NR%8==6 || NR%8==7 || NR%8==0' "$fastq_outputfile".fq > "$fastq_outputfile"_2.fq
