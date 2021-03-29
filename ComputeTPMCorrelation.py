# python3 ComputeTPMCorrelation.py -g unique_ids_10M.pro -t salmon_quant/quant.sf
import readline
import sys, getopt
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def cliParser(argv):
    samgroundinputfile = ''
    decoyinputfile = ''
    proinputfile = ''
    polyesterinputfile = ''
    toolinputfile = ''
    bedinputfile = ''
    fastainputfile1 = ''
    fastainputfile2 = ''
    try:
        opts, args = getopt.getopt(argv, "hs:d:g:p:t:b:1:2:", ["samgroundinputfile", "proinputfile", "polyesterinputfile", "toolinputfile", "bedinputfile", "fastainputfile1", "fastainputfile2"])
    except getopt.GetoptError:
        print('python3 ComputeTPMCorrelation.py -s <samgroundinputfile> -g <proinputfile> -p <polyesterinputfile> -t <toolinputfile> < -b <bedinputfile> / -1 <fastainputfile1> -2 <fastainputfile2> >')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 ComputeTPMCorrelation.py -s <samgroundinputfile> -g <proinputfile> -p <polyesterinputfile> -t <toolinputfile> < -b <bedinputfile> / -1 <fastainputfile1> -2 <fastainputfile2> >')
            sys.exit()
        elif opt in ("-s", "--samgroundinputfile"):
            samgroundinputfile = arg            
        elif opt in ("-g", "--proinputfile"):
            proinputfile = arg
        elif opt in ("-p", "--polyesterinputfile"):
            polyesterinputfile = arg
        elif opt in ("-t", "--toolinputfile"):
            toolinputfile = arg
        elif opt in ("-b", "--bedinputfile"):
            bedinputfile = arg
        elif opt in ["-1", "--fastainputfile1"]:
            fastainputfile1 = arg
        elif opt in ["-2", "--fastainputfile2"]:
            fastainputfile2 = arg
    print ("samgroundinputfile ", samgroundinputfile)
    print ("proinputfile ", proinputfile)
    print ("polyesterinputfile ", polyesterinputfile)
    print ("toolinputfile ", toolinputfile)
    print ("bedinputfile ", bedinputfile)
    print ("fastainputfile1 ", fastainputfile1)
    print ("fastainputfile2 ", fastainputfile2)
    return samgroundinputfile, proinputfile, polyesterinputfile, toolinputfile, bedinputfile, fastainputfile1, fastainputfile2

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
def computeTPMFromBed(proinputfile, bedinputfile):
    print("ground truth from BED.")
    txp_len = dict()
    txp_count = dict()
    with open(proinputfile) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_len[ll[1]] = float(ll[3]) / 1000.0
            txp_count[ll[1]] = 0.0
            line = f.readline()[:-1]

    with open(bedinputfile) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            if ll[0] != "polyA":
                txp_count[ll[3].split(":")[2]] += 1.0
            line = f.readline()[:-1] # read fasta line header
    txp_tpm = computeTPMImpl(list(txp_len.values()), list(txp_count.values()))
    return dict(zip(list(txp_len.keys())[:], txp_tpm[:]))

def computeTPMFromSAM(samgroundinputfile, txp_length):
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
    
    txp_tpm = dict()
    denom = 0.0
    for k in txp_length.keys():
        denom += txp_count[k] / (10**6) / txp_len[k]
    for k in txp_length.keys():
        txp_tpm[k] = txp_count[k] / txp_len[k] / denom
    return txp_tpm

def computeTPMImpl(txp_len, txp_count):
    txp_tpm = []
    denom = 0.0
    for i in range(len(txp_len)):
        denom += txp_count[i] / (10**6) / txp_len[i]
    for i in range(len(txp_len)):
        txp_tpm.append(txp_count[i] / txp_len[i] / denom)
    return txp_tpm

def getTPMFromSalmon(inputfile, isSAMgroundinputfile):
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

def computeCorrelation(groundtruthTPM, toolTPM, isNascent):
    relation = list(set(groundtruthTPM) & set(toolTPM))
    if isNascent:
        relation = list(set(groundtruthTPM) | set(toolTPM))
    print("transcript .pro_count", len(set(groundtruthTPM)), ".salmon count", len(set(toolTPM)))

    name = []
    ground = []
    tool = []
    for k in relation:
        name.append(k)
        if k in groundtruthTPM.keys():
            ground.append(groundtruthTPM[k])
        else:
            ground.append(0.0)
        if k in toolTPM.keys():
            tool.append(toolTPM[k])
        else:
            tool.append(0.0)
    print(name)
    print(ground)
    print(tool)
    df = pd.DataFrame([ground, tool]).transpose()
    df.index = name
    df.columns = ["ground_TPM", "salmon_TPM"]
    print("gene1 30% read count(30%nascent inside), gene2 70% read count(0%nascent inside)")
    print(df)
    print("salmon index transcrpitome", stats.spearmanr(ground, tool))
    print("salmon index transcriptome pearson" , stats.pearsonr(ground, tool))

    # ground = []
    # tool = []
    # for k in relation:
    #     ground.append(groundtruthTPM[k])
    #     tool.append(toolTPM[k])
    # print("transcrpitome+nascent ", stats.spearmanr(ground, tool))
    # print("transcriptome+nascent ", stats.pearsonr(ground, tool))

    # nonascent_ground = []
    # nonascent_tool = []
    # for k in relation:
    #     if k[:4] == "FBtr":
    #         nonascent_ground.append(groundtruthTPM[k])
    #         nonascent_tool.append(toolTPM[k])
    # print("transcriptome ", stats.spearmanr(nonascent_ground, nonascent_tool))
    # print("transcriptome ", stats.pearsonr(nonascent_ground, nonascent_tool))
    
    draw_x = []
    draw_y = []
    for i in range(len(relation)):
        draw_x.append(ground[i])
        draw_y.append(tool[i])
    plt.scatter(draw_x, draw_y)
    a = np.linspace(0.0, 300000.0, num=300001, endpoint=True)
    plt.scatter(a,a, s=1)
    plt.title("x=groundtruthTPM y=toolTPM")
    plt.savefig("TPM.png")

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

if __name__ == "__main__":
    samgroundinputfile, proinputfile, polyesterinputfile, toolinputfile, bedinputfile, fastainputfile1, fastainputfile2 = cliParser(sys.argv[1:])
    isSAMgroundinputfile = False
    if samgroundinputfile:
        isSAMgroundinputfile = True

    toolTPM, txp_length = getTPMFromSalmon(toolinputfile, isSAMgroundinputfile)
        
    isNascent = True
    if samgroundinputfile:
        groundtruthTPM = computeTPMFromSAM(samgroundinputfile, txp_length)
        computeCorrelation(groundtruthTPM, toolTPM, False)
    elif proinputfile:
        if bedinputfile:
            groundtruthTPM = computeTPMFromBed(proinputfile, bedinputfile)
            computeCorrelation(groundtruthTPM, toolTPM, isNascent)
        elif fastainputfile1 and fastainputfile2:
            groundtruthTPM = computeTPMFromFasta(proinputfile, fastainputfile1, fastainputfile2)
            computeCorrelation(groundtruthTPM, toolTPM, isNascent)
        else:
            groundtruthTPM = computeTPMFromPro(proinputfile)
            computeCorrelation(groundtruthTPM, toolTPM, isNascent)
    elif polyesterinputfile:
        groundtruthTPM = computeTPMFromPolyester(polyesterinputfile, txp_length)
        computeCorrelation(groundtruthTPM, toolTPM, isNascent)
    else:
        raise Exception("Neither proinputfile nor polyesterinputfile are given.")
    