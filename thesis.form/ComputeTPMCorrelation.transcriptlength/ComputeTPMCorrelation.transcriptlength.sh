#!/bin/bash

replicate=10
percent=60
basepath="/home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep"
index="jimmy_salmon.txp.add_nascent.multimapping.entropy.1220"
TPMdir="computeTPM"
for i in $(seq 1 ${replicate})
do
    cd ${basepath}${i}/${index}/${index}.idxquant/${TPMdir}/
    python3 /home/scripts/thesis.form/ComputeTPMCorrelation.transcriptlength/ComputeTPMCorrelation.transcriptlength.py -l ../quant.sf -t TPM.txt
    # python3 /home/scripts/thesis.form/ComputeTPM.transcriptlength.py -l ${basepath}${i}/${index}/${index}.idxquant/quant.sf -t ${basepath}${i}/${index}/${index}.idxquant/${TPMdir}/TPM.txt
done