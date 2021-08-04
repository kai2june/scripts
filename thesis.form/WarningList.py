import readline
import sys, getopt
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import operator
import os

def cliParser(argv):
    basepath=''
    nascent_percent=''
    index=''
    try:
        opts, args = getopt.getopt(argv, "hp:n:i:", ["basepath", "nascent_percent", "index"])
    except getopt.GetoptError:
        print('python3 WarningList.py -p <basepath> -n <nascent_percent> -i <index>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 WarningList.py -p <basepath> -n <nascent_percent> -i <index>')
            sys.exit()
        elif opt in ("-p", "--basepath"):
            basepath = arg
        elif opt in ("-n", "--nascent_percent"):
            nascent_percent = arg
        elif opt in ("-i", "--index"):
            index = arg

    # print ("basepath: ", basepath)
    # print ("nascent_percent: ", nascent_percent)
    # print ("index: ", index)

    return basepath, nascent_percent, index

def findFPFNlist(file_prefix, threshold_percent_lower, threshold_percent_upper, least_occurence):
    replicate = 10
    FPFN = dict()
    for i in range(1, replicate+1):
        fp_count = 0 
        with open("/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep" + str(i) + "/salmon.txp/salmon.txp.quant/computeTPM/" + file_prefix + "rank.txt") as f:
            line = f.readline()[:-1]
            while line:
                fp_count += 1
                line = f.readline()[:-1]
        with open("/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep" + str(i) + "/salmon.txp/salmon.txp.quant/computeTPM/" + file_prefix + "rank.txt") as f:
            threshold_lower = math.floor(fp_count * threshold_percent_lower)
            threshold_upper = math.floor(fp_count * threshold_percent_upper) - threshold_lower + 1
            line = f.readline()[:-1]
            while line and threshold_lower > 1:
                line = f.readline()[:-1]
                threshold_lower -= 1
            while line and threshold_upper > 0:
                if not line in FPFN:
                    FPFN[line] = 1
                else:
                    FPFN[line] += 1
                threshold_upper -= 1
                line = f.readline()[:-1]
                
    rlt = []
    for k in FPFN.keys():
        if FPFN[k] >= least_occurence:
            rlt.append(k)
    return rlt
    # intersection = FPFN[1]
    # for k in FPFN.keys():
    #     threshold = math.floor(len(FPFN[k]) * threshold_percent)
    #     intersection = list( set(intersection) & set(FPFN[k][:threshold]) )
    # return intersection


def statEqvClassInSamples(eqv_dir, FPlist, FPratio_threshold):
    txp_id_name = dict()
    stat_eqvclass_occurences = dict()
    stat_eqvclass_FP = dict()
    stat_eqvclass_size = dict()
    for filename in os.listdir(eqv_dir):
        with open(filename) as f:
            line = f.readline()[:-1]
            if not line:
                continue
            txp_count = int(line)
            eqv_count = int(f.readline()[:-1])
            for i in range(txp_count):
                txp_id_name[i] = f.readline()[:-1]
            for i in range(eqv_count):
                txp_list = []
                line = f.readline()[:-1]
                ll = line.split()
                groupSize = int(ll[0])
                FP = 0
                for j in range(1, groupSize+1):
                    txp_name = txp_id_name[int(ll[j])]
                    if txp_name in FPlist:
                        FP += 1
                    txp_list.append(txp_name)
                txp_list = tuple(txp_list)
                if not txp_list in stat_eqvclass_occurences:
                    stat_eqvclass_occurences[txp_list] = 1
                    stat_eqvclass_size[txp_list] = groupSize
                    stat_eqvclass_FP[txp_list] = FP
                else:
                    stat_eqvclass_occurences[txp_list] += 1
                    stat_eqvclass_size[txp_list] += groupSize
                    stat_eqvclass_FP[txp_list] += FP
    ### eqvFPratio
    eqvFPratio = dict()
    for k in stat_eqvclass_occurences.keys():
        eqvFPratio[k] = float(stat_eqvclass_FP[k]) / float(stat_eqvclass_size[k])
    with open("../txp_eqv_FPrank.txt", 'w') as f:
        f.write("transcript false positive: %d\n" %(len(FPlist)))
        for elem in FPlist:
            f.write(str(elem) + "\n")
        f.write("eqvclass false positive ratio more than %f:\n" %FPratio_threshold)
        for counter in range(int(len(eqvFPratio)/2)):
            ky = max(eqvFPratio.items(), key=operator.itemgetter(1))[0]
            if eqvFPratio[ky] == 0.0 or eqvFPratio[ky] < FPratio_threshold:
                break
            for elem in ky:
                f.write(elem + " ")
            f.write(str(float(int(eqvFPratio[ky]*10000))/10000.0) + "\n")
            eqvFPratio.pop(ky, None)

if __name__ == "__main__":
    basepath, nascent_percent, index = cliParser(sys.argv[1:])
    threshold_percent_lower = 0.0
    threshold_percent_upper = 0.5
    least_occurence = 6
    FPintersection = findFPFNlist("FP", threshold_percent_lower, threshold_percent_upper, least_occurence) # 所有sample之中至少在前25%出現一次的有3722個isoform
    print ("threshold_percent_lower, threshold_percent_upper, least_occurence", threshold_percent_lower, threshold_percent_upper, least_occurence)
    print ("size of FPintersection:", len(FPintersection))
    print ("FPintersection", FPintersection)

    FPratio_threshold = 0.3
    statEqvClassInSamples("/home/warningList/eqv/", FPintersection, FPratio_threshold)