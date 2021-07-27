import readline
import sys, getopt
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import operator

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

def findFPFNlist(file_prefix, threshold_percent_lower, threshold_percent_upper, least_occurence, only_report_repi):
    replicate = 10
    FPFN = dict()
    repi = []
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
                if i == only_report_repi:
                    repi.append(line)
                threshold_upper -= 1
                line = f.readline()[:-1]
                
    rlt = []
    for k in FPFN.keys():
        if not k in repi:
            continue
        if FPFN[k] >= least_occurence:
            rlt.append(k)
    return rlt
    # intersection = FPFN[1]
    # for k in FPFN.keys():
    #     threshold = math.floor(len(FPFN[k]) * threshold_percent)
    #     intersection = list( set(intersection) & set(FPFN[k][:threshold]) )
    # return intersection

if __name__ == "__main__":
    basepath, nascent_percent, index = cliParser(sys.argv[1:])
    threshold_percent_lower = 0.0
    threshold_percent_upper = 0.25
    least_occurence = 5
    only_report_repi = 1
    FPintersection = findFPFNlist("FP", threshold_percent_lower, threshold_percent_upper, least_occurence, only_report_repi) # 所有sample之中至少在前25%出現一次的有3722個isoform
    print ("FPintersection", FPintersection)
    print ("size of FPintersection:", len(FPintersection))
    print ("threshold_percent_lower, threshold_percent_upper, least_occurence", threshold_percent_lower, threshold_percent_upper, least_occurence)