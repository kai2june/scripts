#!/bin/bash

replicate=10
# if 0% then 
# {--add_nascent_threshold 1.0 --intron_read_percentage_in_nascent_at_least 1.0
# ./ComputeTPMCorrelation.py -j} => computeTPM.txp/
# elif >0% then 
# {--add_nascent_threshold 0.0 --intron_read_percentage_in_nascent_at_least 0.0
# ./ComputeTPMCorrelation.py -i -j} => computeTPM/
percent=0
tool="jimmy_salmon"
index="txp.add_nascent.multimapping.entropy.1222"
mode="idx"
branch=""
options="-d --transcriptome_size_no_nascent 34611 --add_nascent_threshold 1.0 --intron_read_percentage_in_nascent_at_least 1.0"
# transcript_fasta="gffread.txp.fa"
for i in $(seq 1 ${replicate})
do
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i} && 
    # if [ ! -d jimmy_salmon.${index} ] 
    #     then
    #         mkdir jimmy_salmon.${index}
    # fi
    if [ ! -d jimmy_salmon.${index} ] 
        then
            mkdir jimmy_salmon.${index}
    fi
    cd jimmy_salmon.${index} && 
    if [ ${mode} = "idx" ]
    then
        echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<quant in idx mode"
        # salmon index --genome /home/droso.reference/Drosophila_melanogaster.BDGP6.dna.chromosome.all.fa --input_transcript_gff3 /home/droso.reference/gffread.txp.gff3 --transcripts /home/droso.reference/gffread.txp.fa --output_gene_gff3 jimmy_salmon.gene.gff3 --output_transcript_gene_fasta jimmy_salmon.${index}.fa -i jimmy_salmon.${index}.index -p 60 --keepDuplicates &&
        # salmon quant -i /home/droso.reference/jimmy_salmon.${index}.index -l A -1 ../txpgene${percent}%_1.fasta -2 ../txpgene${percent}%_2.fasta -d -o jimmy_salmon.${index}.${mode}quant --transcriptome_size_no_nascent 34611 --add_nascent_threshold 1.0 --intron_read_percentage_in_nascent_at_least 1.0
        salmon quant -i /home/droso.reference/jimmy_salmon.txp.add_nascent.index -l A -1 ../txpgene${percent}%_1.fasta -2 ../txpgene${percent}%_2.fasta -o jimmy_salmon.${index}.${mode}quant ${options}
    elif [ ${mode} = "aln" ]
    then
        echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<quant in aln mode"
        # salmon quant -t /home/droso.reference/${transcript_fasta} -l A -a ${tool}.${index}.idxquant.sam -d -o ${tool}.${index}.${branch}.${mode}quant
    else
        echo "Exception: neither idx nor aln mode."
        exit -1
    fi
    echo "jimmy_salmon.${index}.${mode}quant" &&
    cd jimmy_salmon.${index}.${mode}quant && 
    rm -rf computeTPM &&
    mkdir computeTPM && 
    cd computeTPM && 
    python3 /home/scripts/ComputeTPMCorrelation.py -g /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%.pro -b /home/0309meeting/0413/txp_vs_genetxp/nascent${percent}%/txpgene${percent}%.rep${i}/txpgene${percent}%.bed -t ../quant.sf -i -j & 
done
