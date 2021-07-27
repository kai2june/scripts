#!/bin/bash

basepath=""
while getopts "p:" argv;
do
    case "${argv}" in
        p)
            basepath=${OPTARG}
            ;;
        *)
            echo "exception: getopts problem."
            ;;
    esac
done

replicate=10
nascent_unique_count=""
all_count=""
for i in $(seq 1 ${replicate})
do
    cd ${basepath}${i}/salmon.txpgene/salmon.txpgene.quant/aux_info
    [ -f eq_classes.txt.gz ] && gunzip eq_classes.txt.gz
    nascent_unique_count=${nascent_unique_count}+`grep '^\<1\>' eq_classes.txt | awk '{if($2 >= 34611) sum+=$NF}END{print sum}'`
    all_count=${all_count}+`cat eq_classes.txt | awk '{sum+=$NF}END{print sum}'`
done

echo "${nascent_unique_count} ${all_count}" | awk -F"+" -v replicate=${replicate} '{for(i=2; i<(NF+1)/2+1; ++i) {j=i+replicate; print $i/$j; sum+=$i/$j;}}END{printf("avg:"sum/replicate"\n")}'