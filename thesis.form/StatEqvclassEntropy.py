import readline
import math
import operator
import numpy as np
from matplotlib import pyplot as plt

def draw_hist(eqvEntropy, category, maxEntropy, y_bound):
    bins = np.arange(start=0, stop=maxEntropy, step=0.1)
    plt.ylim(0,y_bound)
    plt.hist(eqvEntropy, bins=bins, edgecolor='black')
    #plt.xlabel("entropy of weights in an equivalence class", fontsize=16)
    #plt.ylabel("count of the equivalence class", fontsize=16)
    #plt.title("Entropy dist.: Eqv. class with " + category + " #FP/#TP", fontsize=16)
    # plt.title("Entropy distribution: Eqvclass FP/TP ratio " + category + "(" + str(len(eqvEntropy)) + "s)")
    plt.savefig(category + "_entropy" + ".png")
    plt.clf()

def statEqvclassEntropy(eqvfilename, FPfilename, TPTNfilename):
######## record FP, TP transcript
    falsePositive = []
    truePositive = []
    sample_count = 0
    with open(FPfilename) as f:
        sample_count = int(int(f.readline()[:-1]) * 0.5) ### The first line record the count of FP txp.
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
        for i in range(sample_count): ### get TP txp, with the number |FP txp|
            truePositive.append(txp_name)
            line = f.readline()[:-1]
            if not line:
                break
            else:
                txp_name = line.split()[0]
######## record FP/TP ratio and entropy in each equivalence class
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
    return topQuarterFP, bottomQuarterFP, maxEntropy
    # draw_hist(topQuarterFP, "top25%", maxEntropy)
    # draw_hist(bottomQuarterFP, "bottom25%", maxEntropy)

if __name__ == "__main__":
    topQuarter = []
    bottomQuarter = []
    maximumEntropy = 0.0
    for i in range(1, 11):
        eqvfilename = "/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep" + str(i) + "/salmon.txp/salmon.txp.idxquant/aux_info/eq_classes.txt"
        FPfilename = "/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep" + str(i) + "/salmon.txp/salmon.txp.idxquant/computeTPM/FPrank.txt"
        TPTNfilename = "/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep" + str(i) + "/salmon.txp/salmon.txp.idxquant/computeTPM/TPTNrank.txt"
        topQuarterFP, bottomQuarterFP, maxEntropy = statEqvclassEntropy(eqvfilename, FPfilename, TPTNfilename)
        topQuarter.extend(topQuarterFP)
        bottomQuarter.extend(bottomQuarterFP)
        maximumEntropy = max(maxEntropy, maximumEntropy)
    y_bound = 30000
    draw_hist(topQuarter, "higher", maxEntropy, y_bound=y_bound) # higher25%
    draw_hist(bottomQuarter, "lower", maxEntropy, y_bound=y_bound) # lower25%

### @brief 6 minutes execution time for 10 samples about 60000~70000 points plotted
