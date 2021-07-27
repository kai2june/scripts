#!/bin/bash

set -x

cat /home/.bashrc
source /home/.bashrc
cd /mammoth/sailfish_data
/home/sailfish/stage/bin/sailfish quant -i /mammoth/sailfish_data/31 -l IU -1 /mammoth/salmon_data/reads/SRR6269052_1.fastq -2 /mammoth/salmon_data/reads/SRR6269052_2.fastq -o /mammoth/sailfish_data/sailfish_quant_SRR6269052
