import readline
import math
import operator
import numpy as np
from matplotlib import pyplot as plt

def draw_hist(eqvEntropy, category, maxEntropy):
    bins = np.arange(start=0, stop=maxEntropy, step=0.1)
    plt.ylim(0,3000)
    plt.hist(eqvEntropy, bins=bins, edgecolor='black')
    plt.xlabel("entropy in an equivalence class")
    plt.ylabel("count of equivalence class")
    plt.title("Entropy distribution: Eqvclass FP/TP ratio " + category + "(" + str(len(eqvEntropy)) + "s)")
    plt.savefig(category + "_entropy" + ".png")
    plt.clf()

def statEqvclassEntropy(eqvfilename, FPfilename, FNfilename, TPTNfilename):
########
    falsePositive = []
    truePositive = []
    sample_count = 0
    with open(FPfilename) as f:
        sample_count = int(int(f.readline()[:-1]) * 0.5)
        for i in range(sample_count):
            line = f.readline()[:-1]
            falsePositive.append(line.split()[0])
    with open(TPTNfilename) as f:
        line = f.readline()[:-1]
        line = f.readline()[:-1]
        txp_name, txp_tpm = '', 0.0
        while line:
            ll = line.split()
            txp_name, txp_tpm = ll[0], float(ll[1])
            if txp_tpm != 0.0:
                break
            line = f.readline()[:-1]
        for i in range(sample_count):
            truePositive.append(txp_name)
            line = f.readline()[:-1]
            if not line:
                break
            else:
                txp_name = line.split()[0]
########
    txp_id_name = dict()
    eqv_fp_over_tp = dict()
    eqv_entropy = dict()
    txp_count = 0
    eqv_count = 0
    with open(eqvfilename) as f:
        txp_count = int(f.readline()[:-1])
        eqv_count = int(f.readline()[:-1])
        for i in range(txp_count):
            txp_name = f.readline()[:-1]
            txp_id_name[i] = txp_name
        for i in range(eqv_count):
            line = f.readline()[:-1]
            ll = line.split()
            eqv_fp = 0
            eqv_tp = 0
            entropy = 0
            groupSize = int(ll[0])
            for j in range(1, groupSize+1):
                txp_id = int(ll[j])
                if txp_id_name[txp_id] in falsePositive:
                    eqv_fp += 1
                elif txp_id_name[txp_id] in truePositive:
                    eqv_tp += 1
                else:
                    pass
            if eqv_tp == 0:
                if eqv_fp == 0:
                    eqv_fp_over_tp[i] = 0.0
                else:
                    eqv_fp_over_tp[i] = math.inf
            else:
                eqv_fp_over_tp[i] = float(eqv_fp)/float(eqv_tp)
            for j in range(groupSize+1, 2*groupSize+1):
                prob = float(ll[j])
                entropy += -prob*math.log(prob)
            eqv_entropy[i] = entropy
########
    topQuarterFP = []
    bottomQuarterFP = []
    maxEntropy = 0.0
    for k in range(int(eqv_count*0.25)):
        ky = max(eqv_fp_over_tp.items(), key=operator.itemgetter(1))[0]
        topQuarterFP.append(eqv_entropy[ky])
        if eqv_entropy[ky] > maxEntropy:
            maxEntropy = eqv_entropy[ky]
        eqv_fp_over_tp.pop(ky, None)
    for k in range(int(eqv_count*0.5)):
        ky = max(eqv_fp_over_tp.items(), key=operator.itemgetter(1))[0]
        eqv_fp_over_tp.pop(ky, None)
    while len(eqv_fp_over_tp) > 0:
        ky = max(eqv_fp_over_tp.items(), key=operator.itemgetter(1))[0]
        bottomQuarterFP.append(eqv_entropy[ky])
        if eqv_entropy[ky] > maxEntropy:
            maxEntropy = eqv_entropy[ky]
        eqv_fp_over_tp.pop(ky, None)
    draw_hist(topQuarterFP, "top25%", maxEntropy)
    draw_hist(bottomQuarterFP, "bottom25%", maxEntropy)

if __name__ == "__main__":
    eqvfilename = "/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep1/salmon.txp/salmon.txp.idxquant/aux_info/eq_classes.txt"
    FPfilename = "/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep1/salmon.txp/salmon.txp.idxquant/computeTPM/FPrank.txt"
    FNfilename = "/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep1/salmon.txp/salmon.txp.idxquant/computeTPM/FNrank.txt"
    TPTNfilename = "/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep1/salmon.txp/salmon.txp.idxquant/computeTPM/TPTNrank.txt"
    statEqvclassEntropy(eqvfilename, FPfilename, FNfilename, TPTNfilename)