#!/bin/bash

replicate=10
percent=0
tool="jimmy_salmon"
index="txp"
mode="idx"
branch="multimapping_entropy.seqBias.gcBias"
options="-d --seqBias --gcBias"
# options="-d --transcriptome_size_no_nascent 34611 --add_nascent_threshold 1.0 --intron_read_percentage_in_nascent_at_least 1.0"
transcript_fasta="gffread.txp.fa"
for i in $(seq 1 ${replicate})
do
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i} && 
    if [ ! -d ${tool}.${index}.${branch} ] 
        then
            mkdir ${tool}.${index}.${branch}
    fi
    cd ${tool}.${index}.${branch} && 
    if [ ${mode} = "idx" ]
    then
        echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<quant in idx mode"
        salmon quant -i /home/droso.reference/salmon.${index}.index -l A -1 ../txpgene${percent}%_1.fasta -2 ../txpgene${percent}%_2.fasta ${options} -o ${tool}.${index}.${branch}.${mode}quant &&
        echo "dummy line."
    elif [ ${mode} = "aln" ]
    then
        echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<quant in aln mode"
        # salmon quant -t /home/droso.reference/${transcript_fasta} -l A -a ${tool}.${index}.idxquant.sam -d -o ${tool}.${index}.${branch}.${mode}quant
    else
        echo "Exception: neither idx nor aln mode."
        exit -1
    fi
    echo "${tool}.${index}.${branch}.${mode}quant"
    cd ${tool}.${index}.${branch}.${mode}quant && 
    rm -rf computeTPM &&
    mkdir computeTPM && 
    cd computeTPM && 
    python3 /home/scripts/ComputeTPMCorrelation.py -g /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%.pro -b /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%.bed -t ../quant.sf & 
done
