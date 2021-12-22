import readline
import sys, getopt
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import operator

filterout_lower_than = 50.0
basepath = "/home/0309meeting/0413/txp_vs_genetxp/nascent0%/txpgene0%.rep"
basepath2 = "/salmon.txp/salmon.txp.idxquant/computeTPM/TPM.txt"

spearmans = []
for i in range(1,11):
    inputfile = basepath + str(i) + basepath2 
    print(inputfile)   
    ground_arr = []
    estimated_arr = []  
    with open(inputfile, 'r') as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            ground = float(ll[1])
            estimated = float(ll[2])
            if estimated < filterout_lower_than:
                estimated = 0.0
            # print(ground, estimated)
            line = f.readline()[:-1]
            ground_arr.append(ground)
            estimated_arr.append(estimated)

    # print(str(stats.spearmanr(ground_arr, estimated_arr)))
    spearmans.append(float((str(stats.spearmanr(ground_arr, estimated_arr)).split("=")[1].split(",")[0])))
print(spearmans)

sum = 0.0
for elem in spearmans:
    sum += elem
print(sum / len(spearmans))
