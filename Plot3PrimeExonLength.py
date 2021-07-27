import readline
from matplotlib import pyplot as plt
import numpy as np
import re

def PlotintraGene3PrimeExonLengthFPFN(txpGff3, TPM_txt, FPFileName, FNFileName):
    gene_txp = dict()
    txp_gene = dict()
    txp_exonLen = dict()
    txp_tpm = dict()
    txp_FPFN = dict()
    plot_x = []
    plot_y = []
    with open(txpGff3) as f:
        line = f.readline()[:-1]
        while line:
            if line[0] == '#':
                line = f.readline()[:-1]
                continue
            ll = line.split()
            if "geneID" in ll[-1]:
                gene = re.search('geneID=(.+?);', ll[-1]).group(1)
                txp = re.search('ID=(.+?);', ll[-1]).group(1)
                if not gene in gene_txp.keys():
                    gene_txp[gene] = []
                gene_txp[gene].append(txp)
                txp_gene[txp] = gene
            if ll[2] == "exon":
                try:
                    txp = re.search('Parent=(.*)', ll[-1]).group(1)
                    txp_exonLen[txp] = int(ll[4]) - int(ll[3]) + 1
                except AttributeError:
                    txp = ''
            line = f.readline()[:-1]

    with open(TPM_txt) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_tpm[ll[0]] = float(ll[2]) - float(ll[1])
            txp_FPFN[ll[0]] = str()
            line = f.readline()[:-1]
    
    with open(FPFileName) as f:
        line = f.readline()[:-1]
        while line:
            txp_FPFN[line] = "FP"
            line = f.readline()[:-1]
    
    with open(FNFileName) as f:
        line = f.readline()[:-1]
        while line:
            txp_FPFN[line] = "FN"
            line = f.readline()[:-1]
    for k in txp_FPFN.keys():
        if txp_FPFN[k] != "FP" and txp_FPFN[k] != "FN":
            txp_FPFN[k] = "TPTN"

    for k in gene_txp.keys():
        maxtpm_FP = 0.0
        mintpm_FN = 0.0
        FP_exonLength = 0
        FN_exonLength = 0
        TPTN_avgExonLength = 0.0
        count_TPTN = 0.0
        for txp in gene_txp[k]:
            if txp_FPFN[txp] == "FP" and txp_tpm[txp] > maxtpm_FP:
                maxtpm_FP = txp_tpm[txp]
                FP_exonLength = txp_exonLen[txp]
            elif txp_FPFN[txp] == "FN" and txp_tpm[txp] < mintpm_FN:
                mintpm_FN = txp_tpm[txp]
                FN_exonLength = txp_exonLen[txp]
            elif txp_FPFN[txp] == "TPTN":
                TPTN_avgExonLength += txp_exonLen[txp]
                count_TPTN += 1
        if count_TPTN > 0:
            TPTN_avgExonLength /= count_TPTN
        if maxtpm_FP > 0 and mintpm_FN < 0:
            plot_x.append(FP_exonLength)
            plot_y.append(FN_exonLength)
        # if maxtpm_FP > 0 and count_TPTN > 0:
        #     plot_x.append(FP_exonLength)
        #     plot_y.append(TPTN_avgExonLength)
    print(len(plot_x), len(plot_y))

    slope1 = np.arange(4000)
    plt.scatter(plot_x, plot_y, s=3)
    plt.scatter(slope1, slope1, s=2)
    plt.xlabel("FPexonLength")
    plt.ylabel("FNexonLength")
    plt.title("3primeexonLength in each gene")
    plt.savefig("3primeexonLength_in_each_gene.png")

def plot3PrimeExonLength(txpNameTxt, txpGff3):
    txpName = dict()

    with open(txpNameTxt) as f:
        line = f.readline()[:-1]
        while line:
            txpName[line] = 0
            line = f.readline()[:-1]
    with open(txpGff3) as f:
        line = f.readline()[:-1]
        while line:
            if line[0] == '#':
                line = f.readline()[:-1]
                continue
            ll = line.split()
            if ll[2] == "exon":
                name = ll[-1].split('=')[1]
                if name in txpName.keys():
                    txpName[name] = int(ll[4]) - int(ll[3]) + 1
            else:
                line = f.readline()[:-1]
                continue
            line = f.readline()[:-1]
    for k in txpName.keys():
        print(k, txpName[k])

    filename = "3PrimeExonLength"
    bins = np.arange(start=0, stop=2500, step=50)
    val = []
    for k in txpName.keys():
        val.append(txpName[k])
        plt.hist(val, bins=bins, edgecolor='black')
        plt.title(filename + " x-axis:occurence in how many txpGroup, y-axis: txpCount")
        plt.savefig(filename + ".png")
        plt.cla()
        plt.clf()


if __name__ == "__main__":
    # plot3PrimeExonLength("/home/0309meeting/0413/txp_vs_genetxp/txp/0525/salmon/trash.quant/computeTPM/TPTNgroundcount0.txt", 
    #                      "/home/0309meeting/0413/txp_vs_genetxp/txp/Drosophila_melanogaster.BDGP6.80_txp.gff3")

    txpGff3 = "/home/0309meeting/0413/txp_vs_genetxp/txp/Drosophila_melanogaster.BDGP6.80_txp.gff3"
    TPM_txt = "/home/0309meeting/0413/txp_vs_genetxp/txp/0525/salmon/trash.quant/computeTPM/TPM.txt"
    FPFileName = "/home/0309meeting/0413/txp_vs_genetxp/txp/0525/salmon/trash.quant/computeTPM/FPgroundcount0top2001.txt"
    FNFileName = "/home/0309meeting/0413/txp_vs_genetxp/txp/0525/salmon/trash.quant/computeTPM/FNgroundcount0top3229.txt"
    PlotintraGene3PrimeExonLengthFPFN(txpGff3, TPM_txt, FPFileName, FNFileName)