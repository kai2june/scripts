#!/bin/bash

awk 'FNR>0 && NR==FNR {price[$1]=$2; next} FNR>0{printf "%s kallisto%20f truth%20f\n", $1, price[$1], $2*1000000}' /mammoth/kallisto_data/droso1w/Droso_quant/abund.txt /mammoth/flux_simulator_data/droso_1wreads/transcript_fraction.txt > kallisto_ground.txt
