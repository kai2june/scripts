#!/bin/bash  

replicate=10
percent=60
index="txp"
for i in $(seq 1 ${replicate})
do
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/salmon.${index}/salmon.${index}.quant/computeTPM
    mv FPgroundcount* FPrank.txt
    mv FNgroundcount* FNrank.txt
done