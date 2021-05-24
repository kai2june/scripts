# python3 ComputeTPMCorrelation.py -g unique_ids_10M.pro -t salmon_quant/quant.sf
import readline
import sys, getopt
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import operator

def cliParser(argv):
    authorinputfile = ''
    samgroundinputfile = ''
    decoyinputfile = ''
    proinputfile = ''
    polyesterinputfile = ''
    salmoninputfile = ''
    expressinputfile = ''
    bitseqinputfile = ''
    bedinputfile = ''
    fastainputfile1 = ''
    fastainputfile2 = ''
    try:
        opts, args = getopt.getopt(argv, "ha:s:g:p:t:e:q:b:1:2:", ["authorinputfile", "samgroundinputfile", "proinputfile", "polyesterinputfile", "salmoninputfile", "expressinputfile", "bitseqinputfile", "bedinputfile", "fastainputfile1", "fastainputfile2"])
    except getopt.GetoptError:
        print('python3 ComputeTPMCorrelation.py -a <authorinputfile> -s <samgroundinputfile> -g <proinputfile> -p <polyesterinputfile> { -t <salmoninputfile> -e <expressinputfile> -q <bitseqinputfile> } { -b <bedinputfile> / -1 <fastainputfile1> -2 <fastainputfile2> }')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 ComputeTPMCorrelation.py -a <authorinputfile> -s <samgroundinputfile> -g <proinputfile> -p <polyesterinputfile> { -t <salmoninputfile> -e <expressinputfile> -q <bitseqinputfile> } { -b <bedinputfile> / -1 <fastainputfile1> -2 <fastainputfile2> }')
            sys.exit()
        elif opt in ("-a", "--authorinputfile"):
            authorinputfile = arg
        elif opt in ("-s", "--samgroundinputfile"):
            samgroundinputfile = arg            
        elif opt in ("-g", "--proinputfile"):
            proinputfile = arg
        elif opt in ("-p", "--polyesterinputfile"):
            polyesterinputfile = arg
        elif opt in ("-t", "--salmoninputfile"):
            salmoninputfile = arg
        elif opt in ("-e", "--expressinputfile"):
            expressinputfile = arg
        elif opt in ("-q", "--bitseqinputfile"):
            bitseqinputfile = arg
        elif opt in ("-b", "--bedinputfile"):
            bedinputfile = arg
        elif opt in ["-1", "--fastainputfile1"]:
            fastainputfile1 = arg
        elif opt in ["-2", "--fastainputfile2"]:
            fastainputfile2 = arg
    print ("authorinputfile", authorinputfile)
    print ("samgroundinputfile ", samgroundinputfile)
    print ("proinputfile ", proinputfile)
    print ("polyesterinputfile ", polyesterinputfile)
    print ("salmoninputfile ", salmoninputfile)
    print ("expressinputfile ", expressinputfile)
    print ("bitseqinputfile ", bitseqinputfile)
    print ("bedinputfile ", bedinputfile)
    print ("fastainputfile1 ", fastainputfile1)
    print ("fastainputfile2 ", fastainputfile2)
    return authorinputfile, samgroundinputfile, proinputfile, polyesterinputfile, salmoninputfile, expressinputfile, bitseqinputfile, bedinputfile, fastainputfile1, fastainputfile2

def computeTPMFromPro(inputfile):
    print("ground truth from PRO.")
    txp_name = []
    txp_len = []
    txp_count = []
    with open(inputfile) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_name.append(ll[1])
            txp_len.append(float(ll[3]) / 1000.0)
            txp_count.append(float(ll[9]) / 1000.0)
            line = f.readline()[:-1]
    txp_tpm = computeTPMImpl(txp_len, txp_count)
    return dict(zip(txp_name[:], txp_tpm[:]))

def computeTPMFromFasta(proinputfile, fastainputfile1, fastainputfile2):
    print("ground truth from FASTA.")
    txp_len = dict()
    txp_count = dict()
    with open(proinputfile) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_len[ll[1]] = float(ll[3]) / 1000.0
            txp_count[ll[1]] = 0.0
            line = f.readline()[:-1]
    
    with open(fastainputfile1) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split(":")
            txp_count[ll[2]] += 1.0
            line = f.readline()[:-1] # read notused sequence
            line = f.readline()[:-1] # read fasta line header

    with open(fastainputfile2) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split(":")
            txp_count[ll[2]] += 1.0
            line = f.readline()[:-1] # read notused sequence
            line = f.readline()[:-1] # read fasta line header

    txp_tpm = computeTPMImpl(list(txp_len.values()), list(txp_count.values()))
    return dict(zip(list(txp_len.keys())[:], txp_tpm[:]))
def computeTPMFromBed(proinputfile, bedinputfile, isLongerTxpName=False):
    print("ground truth from BED.")
    txp_len = dict()
    txp_count = dict()
    with open(proinputfile) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_name = ll[1] if not isLongerTxpName else ll[1].split(".")[0]
            txp_len[txp_name] = float(ll[3]) / 1000.0
            txp_count[txp_name] = 0.0
            line = f.readline()[:-1]

    with open(bedinputfile) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            if ll[0] != "polyA":
                txp_name = ll[3].split(":")[2] if not isLongerTxpName else ll[3].split(":")[2].split(".")[0]
                txp_count[txp_name] += 1.0
            line = f.readline()[:-1] # read fasta line header
    txp_tpm = computeTPMImpl(list(txp_len.values()), list(txp_count.values()))
    return dict(zip(list(txp_len.keys())[:], txp_tpm[:]))

def computeTPMFromSAM(samgroundinputfile, txp_length):
    print(len(txp_length)) ## sam header has 196354 @SQ  ## 196520 in .gtf
    txp_count = dict()
    for k in txp_length.keys():
        txp_count[k] = 0.0
    with open(samgroundinputfile) as f:
        line = f.readline()[:-1]
        while line[0] == '@':
            line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_count[ll[2]] += 1.0
            line = f.readline()[:-1]
    # print(txp_count)

    txp_tpm = dict()
    denom = 0.0
    for k in txp_length.keys():
        denom += txp_count[k] / (10**6) / txp_length[k]
    for k in txp_length.keys():
        txp_tpm[k] = txp_count[k] / txp_length[k] / denom
    return txp_tpm

def computeTPMImpl(txp_len, txp_count):
    txp_tpm = []
    denom = 0.0
    for i in range(len(txp_len)):
        denom += txp_count[i] / (10**6) / txp_len[i]
    for i in range(len(txp_len)):
        txp_tpm.append(txp_count[i] / txp_len[i] / denom)
    return txp_tpm

def getTPMFromSalmon(inputfile, isSAMgroundinputfile=False):
    txp_tpm = dict()
    txp_length = dict()
    with open(inputfile) as f:
        line = f.readline()[:-1]
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            if isSAMgroundinputfile:
                txp_tpm[ll[0].split(".")[0]] = float(ll[3])
                txp_length[ll[0].split(".")[0]] = float(ll[1]) / 1000.0
            else:
                txp_tpm[ll[0]] = float(ll[3])
                txp_length[ll[0]] = float(ll[1]) / 1000.0
            line = f.readline()[:-1]
    return txp_tpm, txp_length

def getTPMFromExpress(inputfile):
    txp_tpm = dict()
    with open(inputfile) as f:
        line = f.readline()[:-1]
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_tpm[ll[1]] = float(ll[-1])
            line = f.readline()[:-1]
    return txp_tpm
            
def writeTPMFile(groundtruthTPM, toolTPM):
    with open("TPM.txt", 'w') as f:
        for k in toolTPM.keys():
            f.write(k + " " + str(groundtruthTPM[k]) + " " + str(toolTPM[k]) + "\n")

def computeCorrelation(groundtruthTPM, toolTPM, isNascent, convertToLog, proinputfile, onlyComputeProExpressed):
    writeTPMFile(groundtruthTPM, toolTPM)

    relation = list(set(groundtruthTPM) & set(toolTPM))
    if isNascent:
        relation = list(set(groundtruthTPM) | set(toolTPM))
    if onlyComputeProExpressed:
        relation = []
        if proinputfile:
            with open(proinputfile) as f:
                line = f.readline()[:-1]
                while line:
                    ll = line.split()
                    if float(ll[9]) > 0:
                        relation.append(ll[1])
                    line = f.readline()[:-1]

    print("transcript groundtruth_count", len(set(groundtruthTPM)), ".tool count", len(set(toolTPM)))
    print ("{} transcripts account for correlation.".format(len(relation)))            

    decoy = dict()
    if proinputfile:
        with open(proinputfile) as f:
            line = f.readline()[:-1]
            while line:
                if line.split()[1] in relation:
                    decoy[line.split()[1]] = []
                line = f.readline()[:-1]

    col1 = []
    col2 = []
    for k in decoy:
        if k in groundtruthTPM.keys():
            col1.append(groundtruthTPM[k])
        else:
            col1.append(0.0)
        if k in toolTPM.keys():
            col2.append(toolTPM[k])
        else:
            col2.append(0.0)
    if convertToLog:
        for k in decoy.keys():
            if col1 != 0.0:
                col1 = math.log(col1, 2)
            if col2 != 0.0:
                col2 = math.log(col2, 2)
    
    df = pd.DataFrame([col1, col2]).transpose()
    df.index = list(decoy.keys())
    col1_name = ''
    col2_name = ''
    titlename = ''
    figname =''
    if convertToLog:
        col1_name = "ground_logTPM"
        col2_name = "tool_logTPM"
        titlename = "x=groundtruth_logTPM y=tool_logTPM"
        figname = "logTPM.png"
    else:
        col1_name = "ground_TPM"
        col2_name = "tool_TPM"
        titlename = "x=groundtruthTPM y=toolTPM"
        figname = "TPM.png"
    df.columns = [col1_name, col2_name]
    print(df)
    print("salmon index transcrpitome", stats.spearmanr(col1, col2))
    print("salmon index transcriptome pearson" , stats.pearsonr(col1, col2))
    
    draw_x = []
    draw_y = []
    for i in range(len(relation)):
        draw_x.append(col1[i])
        draw_y.append(col2[i])
    plt.scatter(draw_x, draw_y)
    a = np.linspace(0.0, 10000.0, num=10001, endpoint=True)
    plt.scatter(a,a, s=1)
    plt.title(titlename)
    plt.savefig(figname)

    # writeFile of FP, FN, TPTN
    FP = []
    FN = []
    TPTN = []
    FP_dict = dict()
    FN_dict = dict()
    expressed_threshold = 0
    for k in groundtruthTPM.keys():
        if groundtruthTPM[k] <= expressed_threshold and toolTPM[k] > expressed_threshold:
            FP.append(k)
            FP_dict[k] = toolTPM[k]
        elif groundtruthTPM[k] > expressed_threshold and toolTPM[k] <= expressed_threshold:
            FN.append(k)
            FN_dict[k] = toolTPM[k]
        else:
            TPTN.append(k)
    # with open ("FPgroundcount"+str(expressed_threshold)+".txt", 'w') as f:
    #     for k in FP:
    #         f.write(k + '\n')
    # with open ("FNgroundcount"+str(expressed_threshold)+".txt", 'w') as f:
    #     for k in FN:
    #         f.write(k + '\n')
    with open ("TPTNgroundcount"+str(expressed_threshold)+".txt", 'w') as f:
        for k in TPTN:
            f.write(k + '\n')

    # get topN FP
    topN_FP = len(FP_dict)
    FPtopN = []
    for counter in range(topN_FP):
        ky = max(FP_dict.items(), key=operator.itemgetter(1))[0]
        FPtopN.append(ky)
        FP_dict.pop(ky, None)
    with open("FPgroundcount"+str(expressed_threshold)+"top"+str(topN_FP)+".txt", 'w') as f:
        for k in FPtopN:
            f.write(k + '\n')

    topN_FN = len(FN_dict)
    FNtopN = []
    for counter in range(topN_FN):
        ky = max(FN_dict.items(), key=operator.itemgetter(1))[0]
        FNtopN.append(ky)
        FN_dict.pop(ky, None)
    with open("FNgroundcount"+str(expressed_threshold)+"top"+str(topN_FN)+".txt", 'w') as f:
        for k in FNtopN:
            f.write(k + '\n')

def computeTPMFromPolyester(polyesterinputfile, txp_length):
    txp_tpm = dict()
    denom = 0.0
    with open(polyesterinputfile) as f:
        line = f.readline()[:-1]
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            ll2 = ll[0].split("\"")
            print(ll, ll2[1])
            txp_tpm[ll2[1]] = float(ll[-1]) / float(txp_length[ll2[1]])
            denom += float(ll[-1])*(10**-6) / float(txp_length[ll2[1]])
            line = f.readline()[:-1]
    for k in txp_tpm.keys():
        txp_tpm[k] /= denom

    return txp_tpm

def getTPMFromAuthor(authorinputfile):
    txp_tpm = dict()
    with open(authorinputfile) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_tpm[ll[0]] = float(ll[1])
            line = f.readline()[:-1]
    return txp_tpm

def getRPKMFromBitseq(bitseqinputfile):
    txp_rpkm = dict()
    with open(bitseqinputfile) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_rpkm[ll[0]] = float(ll[1])
            line = f.readline()[:-1]
    return txp_rpkm

if __name__ == "__main__":
    authorinputfile, samgroundinputfile, proinputfile, polyesterinputfile, salmoninputfile, expressinputfile, bitseqinputfile, bedinputfile, fastainputfile1, fastainputfile2 = cliParser(sys.argv[1:])
    isSAMgroundinputfile = False
    if samgroundinputfile:
        isSAMgroundinputfile = True
    
    salmonTPM = ''
    txp_length = ''
    expressTPM = ''
    bitseqRPKM = ''
    if salmoninputfile:
        salmonTPM, txp_length = getTPMFromSalmon(salmoninputfile, isSAMgroundinputfile)
    elif expressinputfile:
        expressTPM = getTPMFromExpress(expressinputfile)
    elif bitseqinputfile:
        bitseqRPKM = getRPKMFromBitseq(bitseqinputfile)
    
    isNascent = True
    isLongerTxpName = True
    if authorinputfile and bitseqinputfile:
        groundtruthTPM = getTPMFromAuthor(authorinputfile)
        computeCorrelation(groundtruthTPM, salmonTPM, False, True)
    elif proinputfile and bedinputfile and bitseqinputfile:
        groundtruthTPM = computeTPMFromBed(proinputfile, bedinputfile, isLongerTxpName=True)
        computeCorrelation(groundtruthTPM, bitseqRPKM, False, convertToLog=True)
    elif samgroundinputfile:
        groundtruthTPM = computeTPMFromSAM(samgroundinputfile, txp_length)
        computeCorrelation(groundtruthTPM, salmonTPM, False, True)
    elif proinputfile and bedinputfile:
        groundtruthTPM = computeTPMFromBed(proinputfile, bedinputfile, isLongerTxpName=True)
        if salmoninputfile:
            computeCorrelation(groundtruthTPM, salmonTPM, isNascent=True, convertToLog=False, proinputfile=proinputfile, onlyComputeProExpressed=False)
        elif expressinputfile:
            computeCorrelation(groundtruthTPM, expressTPM, isNascent=True, convertToLog=False, proinputfile=proinputfile, onlyComputeProExpressed=False)
    elif proinputfile and fastainputfile1 and fastainputfile2:
        groundtruthTPM = computeTPMFromFasta(proinputfile, fastainputfile1, fastainputfile2)
        computeCorrelation(groundtruthTPM, salmonTPM)
    elif proinputfile:
        groundtruthTPM = computeTPMFromPro(proinputfile)
        computeCorrelation(groundtruthTPM, salmonTPM)
    elif polyesterinputfile:
        groundtruthTPM = computeTPMFromPolyester(polyesterinputfile, txp_length)
        computeCorrelation(groundtruthTPM, salmonTPM)
    else:
        raise Exception("Neither proinputfile nor polyesterinputfile are given.")
    
