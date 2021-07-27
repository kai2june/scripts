import sys, getopt
import readline
import matplotlib.pyplot as plt
import numpy as np
import operator

def cliParser(argv):
    eq_classes_txt = ''
    FP_txt = ''
    FN_txt = ''
    TPTN_txt = ''
    try:
        opts, args = getopt.getopt(argv, "he:p:n:t:", ["eq_classes_txt", "FP_txt", "FN_txt", "TPTN_txt"])
    except getopt.GetoptError:
        print('python3 ComputeTPMCorrelation.py -e <eq_classes_txt> -p <FP_txt> -n <FN_txt> -t <TPTN_txt>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 ComputeTPMCorrelation.py -e <eq_classes_txt> -p <FP_txt> -n <FN_txt> -t <TPTN_txt>')
            sys.exit()
        elif opt in ("-e", "--eq_classes_txt"):
            eq_classes_txt = arg
        elif opt in ("-p", "--FP_txt"):
            FP_txt = arg
        elif opt in ("-n", "--FN_txt"):
            FN_txt = arg
        elif opt in ("-t", "--TPTN_txt"):
            TPTN_txt = arg

    print ("eq_classes_txt: ", eq_classes_txt)
    print ("FP_txt: ", FP_txt)
    print ("FN_txt: ", FN_txt)
    print ("TPTN_txt: ", TPTN_txt)

    return eq_classes_txt, FP_txt, FN_txt, TPTN_txt

def statOccurence(eq_classes_txt, FP_txt, FN_txt, TPTN_txt):
    with open(eq_classes_txt) as f:
        txp_count = int(f.readline()[:-1])
        txpGroup_count = int(f.readline()[:-1])
        txpid_txpname = dict()
        txpname_occInTxpGroup = dict()

        id = 0
        for i in range(txp_count):
            txpid_txpname[id] = f.readline()[:-1]
            txpname_occInTxpGroup[txpid_txpname[id]] = 0
            id += 1
        for i in range(txpGroup_count):
            line = f.readline()[:-1]
            ll = line.split()
            for j in range(1, int(ll[0])+1):
                txpname_occInTxpGroup[txpid_txpname[int(ll[j])]] += 1

    # sort FP_dict
    FP_dict = dict()
    with open(FP_txt) as f:
        line = f.readline()[:-1]
        while line:
            FP_dict[line] = txpname_occInTxpGroup[line]
            line = f.readline()[:-1]
    topN = len(FP_dict)
    FPtopN = []
    FPtopN_occ = [] 
    for counter in range(topN):
        ky = max(FP_dict.items(), key=operator.itemgetter(1))
        FPtopN.append(ky[0])
        FPtopN_occ.append(ky[1])
        FP_dict.pop(ky[0], None)
    with open("FPtop"+str(topN)+".txpGroup.txt", 'w') as f:
        for i in range(len(FPtopN)):
            f.write(str(FPtopN[i]) + " " + str(FPtopN_occ[i]) + '\n')

    def draw_hist(filename):
        bins = np.arange(start=0, stop=30, step=5)
        with open(filename) as f:
            val = []
            line = f.readline()[:-1]
            while line:
                val.append(txpname_occInTxpGroup[line])
                line = f.readline()[:-1]
            plt.hist(val, bins=bins, edgecolor='black')
            plt.title(filename + " x-axis:occurence in how many txpGroup, y-axis: txpCount")
            plt.savefig(filename + ".png")
            plt.cla()
            plt.clf()
    draw_hist(FP_txt)
    # draw_hist(FN_txt)
    # draw_hist(TPTN_txt)

if __name__ == "__main__":
    eq_classes_txt, FP_txt, FN_txt, TPTN_txt = cliParser(sys.argv[1:])
    statOccurence(eq_classes_txt, FP_txt, FN_txt, TPTN_txt)