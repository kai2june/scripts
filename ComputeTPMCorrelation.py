# python3 ComputeTPMCorrelation.py -g unique_ids_10M.pro -t salmon_quant/quant.sf
import readline
import sys, getopt
from scipy import stats

def cliParser(argv):
    groundtruthinputfile = ''
    toolinputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hg:t:", ["groundtruthinputfile", "toolinputfile"])
    except getopt.GetoptError:
        print('python3 ComputeTPMCorrelation.py -g <groundtruthinputfile> -t <toolinputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 ComputeTPMCorrelation.py -g <groundtruthinputfile> -t <toolinputfile>')
            sys.exit()
        elif opt in ("-g", "--groundtruthinputfile"):
            groundtruthinputfile = arg
        elif opt in ("-t", "--toolinputfile"):
            toolinputfile = arg
    print ("groundtruthinputfile ", groundtruthinputfile)
    print ("toolinputfile ", toolinputfile)
    return groundtruthinputfile, toolinputfile

def computeTPMFromPro(inputfile):
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

def computeTPMImpl(txp_len, txp_count):
    txp_tpm = []
    denom = 0.0
    for i in range(len(txp_len)):
        denom += txp_count[i] / (10**6) / txp_len[i]
    for i in range(len(txp_len)):
        txp_tpm.append(txp_count[i] / txp_len[i] / denom)
    return txp_tpm

def getTPMFromTxt(inputfile):
    txp_tpm = dict()
    with open(inputfile) as f:
        line = f.readline()[:-1]
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_tpm[ll[0]] = float(ll[3])
            line = f.readline()[:-1]
    return txp_tpm

def computeCorrelation(groundtruthTPM, toolTPM):
    intersect = list(set(groundtruthTPM) & set(toolTPM))
    
    ground = []
    tool = []
    for k in intersect:
        ground.append(groundtruthTPM[k])
        tool.append(toolTPM[k])
    print("transcrpitome+nascent ", stats.spearmanr(ground, tool))

    nonascent_ground = []
    nonascent_tool = []
    for k in intersect:
        if k[:4] == "FBtr":
            nonascent_ground.append(groundtruthTPM[k])
            nonascent_tool.append(toolTPM[k])
    print("transcriptome ", stats.spearmanr(nonascent_ground, nonascent_tool))

if __name__ == "__main__":
    groundtruthinputfile, toolinputfile = cliParser(sys.argv[1:])
    groundtruthTPM = computeTPMFromPro(groundtruthinputfile)
    toolTPM = getTPMFromTxt(toolinputfile)
    computeCorrelation(groundtruthTPM, toolTPM)
