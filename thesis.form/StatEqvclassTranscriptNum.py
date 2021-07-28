import readline
import numpy as np
from matplotlib import pyplot as plt

def draw_hist(eqvCount, category, maxGroupSize):
    bins = np.arange(start=0, stop=maxGroupSize, step=1)
    plt.ylim(0, 7000)
    plt.hist(eqvCount, bins=bins, edgecolor='black')
    plt.xlabel("average transcript number in an equivalence class")
    plt.ylabel("count of transcript")
    if category == "FalsePositive":
        plt.title("Eqvclass size distribution: " + category + " top50%(" + str(len(eqvCount)) + "s)")
    else:
        plt.title("Eqvclass size distritution: " + category + " " + str(len(eqvCount)) + "s")
    plt.savefig(category + "_avgEqvTranscriptNum" + ".png")
    plt.clf()

def statEqvclasstTranscriptNum(eqvfilename, FPfilename, FNfilename, TPTNfilename):
    txp_id_name = dict()
    txp_name_eqvTranscriptNum = dict()
    txp_name_eqvCount = dict()
    maxGroupSize = 0
    with open(eqvfilename) as f:
        txp_count = int(f.readline()[:-1])
        eqv_count = int(f.readline()[:-1])
        for i in range(txp_count):
            txp_name = f.readline()[:-1]
            txp_id_name[i] = txp_name
            txp_name_eqvTranscriptNum[txp_name] = 0
            txp_name_eqvCount[txp_name] = 0

        for i in range(eqv_count):
            line = f.readline()[:-1]
            ll = line.split()
            for j in range(1, int(ll[0])+1):
                if int(ll[0]) > maxGroupSize:
                    maxGroupSize = int(ll[0])
                txp_name = txp_id_name[int(ll[j])]
                txp_name_eqvTranscriptNum[txp_name] += int(ll[0])
                txp_name_eqvCount[txp_name] += 1
    
    sample_count = 0
    FP_eqvCount = []
    with open(FPfilename) as f:
        sample_count = int(int(f.readline()[:-1]) * 0.5)
        for i in range(sample_count):
            line = f.readline()[:-1]
            if not line:
                break
            else:
                txp_name = line.split()[0]
            FP_eqvCount.append(float(txp_name_eqvTranscriptNum[txp_name]) / float(txp_name_eqvCount[txp_name]))

    TP_eqvCount = []
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
            TP_eqvCount.append(float(txp_name_eqvTranscriptNum[txp_name]) / float(txp_name_eqvCount[txp_name]))
            line = f.readline()[:-1]
            if not line:
                break
            else:
                txp_name = line.split()[0]
    return FP_eqvCount, TP_eqvCount, maxGroupSize

if __name__ == "__main__":
    FP_total = []
    TP_total = []
    sz = 0
    for i in range(1, 11):
        eqvfilename = "/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep" + str(i) + "/salmon.txp/salmon.txp.idxquant/aux_info/eq_classes.txt"
        FPfilename = "/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep" + str(i) + "/salmon.txp/salmon.txp.idxquant/computeTPM/FPrank.txt"
        FNfilename = "/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep" + str(i) + "/salmon.txp/salmon.txp.idxquant/computeTPM/FNrank.txt"
        TPTNfilename = "/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep" + str(i) + "/salmon.txp/salmon.txp.idxquant/computeTPM/TPTNrank.txt"
        FP_eqvCount, TP_eqvCount, maxGroupSize = statEqvclasstTranscriptNum(eqvfilename, FPfilename, FNfilename, TPTNfilename)
        FP_total = FP_total + FP_eqvCount
        TP_total = TP_total + TP_eqvCount
        if maxGroupSize > sz:
            sz = maxGroupSize
    draw_hist(FP_total, "FalsePositive", sz)
    draw_hist(TP_total, "TruePositive", sz)
