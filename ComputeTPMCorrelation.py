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
    kallistoinputfile = ''
    cufflinksinputfile = ''
    expressinputfile = ''
    bitseqinputfile = ''
    rseminputfile = ''
    bedinputfile = ''
    fastainputfile1 = ''
    fastainputfile2 = ''
    falseNegativeinputfile = ''
    isUnion=False # False:only transcript are used for correlation, True:transcripts+others are used for correlation
    isNascent=False # False: nascent aren't used for correlation, True:nascent are used for correlation
    try:
        opts, args = getopt.getopt(argv, "ha:s:g:p:t:k:c:e:q:m:b:1:2:n:ij", ["authorinputfile", "samgroundinputfile", "proinputfile", "polyesterinputfile", "salmoninputfile", "kallistoinputfile", "cufflinksinputfile", "expressinputfile", "bitseqinputfile", "rseminputfile", "bedinputfile", "fastainputfile1", "fastainputfile2", "falseNegativeinputfile", "isUnion", "isNascent"])
    except getopt.GetoptError:
        print('python3 ComputeTPMCorrelation.py -a <authorinputfile> -s <samgroundinputfile> -g <proinputfile> -p <polyesterinputfile> { -t <salmoninputfile> -k <kallistoinputfile> -c <cufflinksinputfile> -e <expressinputfile> -q <bitseqinputfile> -m <rseminputfile>} { -b <bedinputfile> / -1 <fastainputfile1> -2 <fastainputfile2> } -n <falseNegativeinputfile> -i -j')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 ComputeTPMCorrelation.py -a <authorinputfile> -s <samgroundinputfile> -g <proinputfile> -p <polyesterinputfile> { -t <salmoninputfile> -k <kallistoinputfile> -c <cufflinksinputfile> -e <expressinputfile> -q <bitseqinputfile> -m <rseminputfile>} { -b <bedinputfile> / -1 <fastainputfile1> -2 <fastainputfile2> -n <falseNegativeinputfile>} -i -j')
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
        elif opt in ("-k", "--kallistoinputfile"):
            kallistoinputfile = arg
        elif opt in ("-c", "--cufflinksinputfile"):
            cufflinksinputfile = arg
        elif opt in ("-e", "--expressinputfile"):
            expressinputfile = arg
        elif opt in ("-q", "--bitseqinputfile"):
            bitseqinputfile = arg
        elif opt in ("-m", "--rseminputfile"):
            rseminputfile = arg
        elif opt in ("-b", "--bedinputfile"):
            bedinputfile = arg
        elif opt in ["-1", "--fastainputfile1"]:
            fastainputfile1 = arg
        elif opt in ["-2", "--fastainputfile2"]:
            fastainputfile2 = arg
        elif opt in ["-n", "--falseNegativeinputfile"]:
            falseNegativeinputfile = arg
        elif opt in ["-i", "--isUnion"]:
            isUnion = True
        elif opt in ["-j", "--isNascent"]:
            isNascent = True
    print ("authorinputfile", authorinputfile)
    print ("samgroundinputfile ", samgroundinputfile)
    print ("proinputfile ", proinputfile)
    print ("polyesterinputfile ", polyesterinputfile)
    print ("salmoninputfile/sailfishinputfile ", salmoninputfile)
    print ("kallistoinputfile ", kallistoinputfile)
    print ("cufflinksinputfile ", cufflinksinputfile)
    print ("expressinputfile ", expressinputfile)
    print ("bitseqinputfile ", bitseqinputfile)
    print ("rseminputfile", rseminputfile)
    print ("bedinputfile ", bedinputfile)
    print ("fastainputfile1 ", fastainputfile1)
    print ("fastainputfile2 ", fastainputfile2)
    print ("falseNegativeinputfile ", falseNegativeinputfile)
    print ("isUnion ", isUnion)
    print ("isNascent ", isNascent)
    return authorinputfile, samgroundinputfile, proinputfile, polyesterinputfile, salmoninputfile, kallistoinputfile, cufflinksinputfile, expressinputfile, bitseqinputfile, rseminputfile, bedinputfile, fastainputfile1, fastainputfile2, falseNegativeinputfile, isUnion, isNascent

def writeTPMFile(groundtruthTPM, toolTPM):
    with open("TPM.txt", 'w') as f:
        if len(toolTPM.keys()) > len(groundtruthTPM.keys()):
            for k in toolTPM.keys():
                if not k in groundtruthTPM.keys():
                    f.write(k + " " + str(0.0) + " " + str(toolTPM[k]) + " NotInGroundFile" + "\n")
                else:
                    f.write(k + " " + str(groundtruthTPM[k]) + " " + str(toolTPM[k]) + "\n")
        else:
            for k in groundtruthTPM.keys():
                if not k in toolTPM.keys():
                    f.write(k + " " + str(groundtruthTPM[k]) + " " + str(0.0) + " NotInToolIndex" + "\n")
                else:
                    f.write(k + " " + str(groundtruthTPM[k]) + " " + str(toolTPM[k]) + "\n" )

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

### ground truth from bed rather than pro to account for not including polyA reads 
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

def getTPMFromKallisto(inputfile):
    txp_tpm = dict()
    with open(inputfile) as f:
        line = f.readline()[:-1]
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_tpm[ll[0]] = float(ll[-1])
            line = f.readline()[:-1]
    return txp_tpm

def getTPMFromExpress(inputfile):
    txp_tpm = dict()
    with open(inputfile) as f:
        line = f.readline()[:-1]
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_tpm[ll[1].split('.')[0]] = float(ll[-1])
            line = f.readline()[:-1]
    return txp_tpm

def getFPKMFromCufflinksThenConvertToTPM(inputfile):
    txp_tpm = dict()
    sum = 0.0
    with open(inputfile) as f:
        line = f.readline()[:-1]
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_tpm[ll[0]] = float(ll[9])
            sum += float(ll[9])
            line = f.readline()[:-1]
    for k in txp_tpm.keys():
        txp_tpm[k] = txp_tpm[k] * (10**6) / sum
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

def getTPMFromrsem(rseminputfile):
    txp_rpkm = dict()
    with open(rseminputfile) as f:
        line = f.readline()[:-1]
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_rpkm[ll[0]] = float(ll[5])
            line = f.readline()[:-1]
    return txp_rpkm

def writeFalseNegativeList(groundtruthTPM, toolTPM, filename):
    FN = []
    for k in groundtruthTPM.keys():
        if groundtruthTPM[k] > 0 and toolTPM[k] == 0:
            FN.append(k)
    with open(filename, 'w') as f:
        for k in FN:
            f.write(k + '\n')
    return FN

def readFalseNegativeList(filename):
    FN = []
    with open(filename) as f:
        line = f.readline()[:-1]
        while line:
            FN.append(line)
            line = f.readline()[:-1]
    return FN

def writeRank(groundtruthTPM, toolTPM):
    print ("Writing FP,FN,TPTN rank... .txt")
    # writeFile of FP, FN, TPTN
    relation = list(set(groundtruthTPM) & set(toolTPM))
    FP_dict = dict()
    FN_dict = dict()
    TPTN_dict = dict()
    expressed_threshold = 0
    for k in relation:
        if groundtruthTPM[k] <= expressed_threshold and toolTPM[k] > expressed_threshold:
            FP_dict[k] = toolTPM[k] - 0
        elif groundtruthTPM[k] > expressed_threshold and toolTPM[k] <= expressed_threshold:
            FN_dict[k] = 0 - groundtruthTPM[k]
        else:
            TPTN_dict[k] = abs(groundtruthTPM[k] - toolTPM[k])

    # get topN FP
    with open("FPrank.txt", 'w') as f:
        f.write(str(len(FP_dict)) + '\n')
        for counter in range(len(FP_dict)):
            ky = max(FP_dict.items(), key=operator.itemgetter(1))[0]
            f.write(ky + " " + str(FP_dict[ky]) + '\n')
            FP_dict.pop(ky, None)

    with open("FNrank.txt", 'w') as f:
        f.write(str(len(FN_dict)) + '\n')
        for counter in range(len(FN_dict)):
            ky = min(FN_dict.items(), key=operator.itemgetter(1))[0]
            f.write(ky + " " + str(FN_dict[ky]) + '\n')
            FN_dict.pop(ky, None)

    with open("TPTNrank.txt", 'w') as f:
        f.write(str(len(TPTN_dict)) + '\n')
        for counter in range(len(TPTN_dict)):
            ky = min(TPTN_dict.items(), key=operator.itemgetter(1))[0]
            f.write(ky + " " + str(TPTN_dict[ky]) + '\n')
            TPTN_dict.pop(ky, None)

def computeCorrelation(groundtruthTPM, toolTPM, isUnion, convertToLog, proinputfile, onlyComputeProExpressed, removeTxpList):
    writeTPMFile(groundtruthTPM, toolTPM)

    relation = list(set(groundtruthTPM) & set(toolTPM))
    if isUnion:
        relation = list(set(groundtruthTPM) | set(toolTPM))
    if onlyComputeProExpressed and proinputfile:
        relation = []
        with open(proinputfile) as f:
            line = f.readline()[:-1]
            while line:
                ll = line.split()
                if float(ll[9]) > 0:
                    relation.append(ll[1])
                line = f.readline()[:-1]
    if removeTxpList:
        for k in removeTxpList:
            if k in relation:
                relation.remove(k)

    # decoy = dict()
    # if proinputfile:
    #     with open(proinputfile) as f:
    #         line = f.readline()[:-1]
    #         while line:
    #             if line.split()[1] in relation:
    #                 decoy[line.split()[1]] = []
    #             line = f.readline()[:-1]

    print("transcript groundtruth_count", len(set(groundtruthTPM)), ".tool count", len(set(toolTPM)))
    if not isUnion:
        print ("{} intersection transcripts account for correlation.".format(len(relation)))         
    else:
        print ("{} union transcripts account for correlation.".format(len(relation)))         


    col1 = []
    col2 = []
    for k in relation:
        if k in groundtruthTPM.keys():
            col1.append(groundtruthTPM[k])
        else:
            col1.append(0.0)
        if k in toolTPM.keys():
            col2.append(toolTPM[k])
        else:
            col2.append(0.0)
    if convertToLog:
        for i in range(len(col1)):
            if col1[i] != 0.0:
                col1[i] = math.log(col1[i], 2)
            if col2[i] != 0.0:
                col2[i] = math.log(col2[i], 2)
    
    df = pd.DataFrame([col1, col2]).transpose()
    df.index = relation
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

    with open("correlation.txt", 'w') as f:
        spearman = str(stats.spearmanr(col1, col2)).split("=")[1].split(",")[0]
        pearson = str(stats.pearsonr(col1, col2)).split("(")[1].split(",")[0]
        f.write("spearman {}\n".format(spearman))
        f.write("pearson {}\n".format(pearson))

    writeRank(groundtruthTPM, toolTPM)
    
    # print(len(col1), len(col2))
    # draw_x = []
    # draw_y = []
    # for i in range(len(col1)):
    #     draw_x.append(col1[i])
    #     draw_y.append(col2[i])
    # plt.scatter(draw_x, draw_y)
    # a = np.linspace(0.0, 10000.0, num=10001, endpoint=True)
    # plt.scatter(a,a, s=1)
    # plt.title(titlename)
    # plt.savefig(figname)

if __name__ == "__main__":
    authorinputfile, samgroundinputfile, proinputfile, polyesterinputfile, salmoninputfile, kallistoinputfile, cufflinksinputfile, expressinputfile, bitseqinputfile, rseminputfile, bedinputfile, fastainputfile1, fastainputfile2, falseNegativeinputfile, isUnion, isNascent = cliParser(sys.argv[1:])

    isSAMgroundinputfile = False
    if samgroundinputfile:
        isSAMgroundinputfile = True
    
    salmonTPM = ''
    txp_length = ''
    kallistoTPM = ''
    cufflinksTPM = ''
    expressTPM = ''
    bitseqRPKM = ''
    rsemTPM = ''
    if salmoninputfile:
        salmonTPM, txp_length = getTPMFromSalmon(salmoninputfile, isSAMgroundinputfile)
    elif kallistoinputfile:
        kallistoTPM = getTPMFromKallisto(kallistoinputfile)
    elif cufflinksinputfile:
        cufflinksTPM = getFPKMFromCufflinksThenConvertToTPM(cufflinksinputfile)
    elif expressinputfile:
        expressTPM = getTPMFromExpress(expressinputfile)
    elif bitseqinputfile:
        bitseqRPKM = getRPKMFromBitseq(bitseqinputfile)
    elif rseminputfile:
        rsemTPM = getTPMFromrsem(rseminputfile)

    convertToLog = False
    onlyComputeProExpressed = False
    isLongerTxpName = True
    removeTxpList = []
    # falseNegativeinputfile = '/home/0309meeting/0413/txp_vs_genetxp/txp/0601/aligner.bowtie2.quant/salmonFN.txt'
    # falseNegativeinputfile = '/home/0309meeting/0413/txp_vs_genetxp/txp/0601/RSEM.bowtie2/rsem.quant.stat/rsemFN.txt'

    if authorinputfile:
        groundtruthTPM = getTPMFromAuthor(authorinputfile)
        if expressinputfile:
            computeCorrelation(groundtruthTPM, expressTPM, isUnion, convertToLog, proinputfile, onlyComputeProExpressed, removeTxpList)
        elif bitseqinputfile:
            computeCorrelation(groundtruthTPM, bitseqRPKM, False, True)
    elif proinputfile and bedinputfile and bitseqinputfile:
        groundtruthTPM = computeTPMFromBed(proinputfile, bedinputfile, isLongerTxpName=True)
        computeCorrelation(groundtruthTPM, bitseqRPKM, False, convertToLog=True)
    elif samgroundinputfile:
        groundtruthTPM = computeTPMFromSAM(samgroundinputfile, txp_length)
        computeCorrelation(groundtruthTPM, salmonTPM, False, True)
    elif proinputfile and bedinputfile:
        groundtruthTPM = computeTPMFromBed(proinputfile, bedinputfile, isLongerTxpName=True)
        if salmoninputfile:
            if falseNegativeinputfile:
                writeFalseNegativeList(groundtruthTPM, salmonTPM, "salmonFN.txt")
                removeTxpList = readFalseNegativeList(falseNegativeinputfile)
            if not isNascent:
                for k in groundtruthTPM.keys():
                    if "FBgn" in k:
                        removeTxpList.append(k)
            computeCorrelation(groundtruthTPM, salmonTPM, isUnion=isUnion, convertToLog=convertToLog, proinputfile=proinputfile, onlyComputeProExpressed=onlyComputeProExpressed, removeTxpList=removeTxpList)
        elif expressinputfile:
            computeCorrelation(groundtruthTPM, expressTPM, isUnion=isUnion, convertToLog=convertToLog, proinputfile=proinputfile, onlyComputeProExpressed=onlyComputeProExpressed, removeTxpList=removeTxpList)
        elif rseminputfile:
            if falseNegativeinputfile:
                writeFalseNegativeList(groundtruthTPM, rsemTPM, "rsemFN.txt")
                removeTxpList = readFalseNegativeList(falseNegativeinputfile)
            computeCorrelation(groundtruthTPM, rsemTPM, isUnion=isUnion, convertToLog=convertToLog, proinputfile=proinputfile, onlyComputeProExpressed=onlyComputeProExpressed, removeTxpList=removeTxpList)
        elif kallistoinputfile:
            computeCorrelation(groundtruthTPM, kallistoTPM, isUnion=isUnion, convertToLog=convertToLog, proinputfile=proinputfile, onlyComputeProExpressed=onlyComputeProExpressed, removeTxpList=removeTxpList)
        elif cufflinksinputfile:
            computeCorrelation(groundtruthTPM, cufflinksTPM, isUnion=isUnion, convertToLog=convertToLog, proinputfile=proinputfile, onlyComputeProExpressed=onlyComputeProExpressed, removeTxpList=removeTxpList)
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
    
