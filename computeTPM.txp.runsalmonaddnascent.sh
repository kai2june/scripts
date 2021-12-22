replicate=10
for i in $(seq 1 ${replicate})
do
    cd /home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep${i}/jimmy_salmon.txp.add_nascent.multimapping.entropy.1222/jimmy_salmon.txp.add_nascent.multimapping.entropy.1222.idxquant
    mkdir computeTPM.txp 
    cd computeTPM.txp
    python3 /home/scripts/ComputeTPMCorrelation.py -g /home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep${i}/txpgene0%.pro -b /home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep${i}/txpgene0%.bed -t ../quant.sf -j & 
done