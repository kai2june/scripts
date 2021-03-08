#!/bin/bash
for i in $(seq 49 52);
do
	prefetch SRR62690${i}
        fastq-dump --split-files SRR62690${i}
done	
