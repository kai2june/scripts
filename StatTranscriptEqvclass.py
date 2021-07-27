import readline
import sys, getopt
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import operator

def cliParser(argv):
    eqvclassinputfile = ''
    try:
        opts, args = getopt.getopt(argv, "he:", ["eqvclassinputfile"])
    except getopt.GetoptError:
        print('python3 StatTranscriptEqvclass.py -e <eqvclassinputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 StatTranscriptEqvclass.py -e <eqvclassinputfile>')
            sys.exit(3)
        elif opt in ("-e", "--eqvclassinputfile"):
            eqvclassinputfile = arg
    print("eqvclassinputfile ", eqvclassinputfile)
    return eqvclassinputfile

def eqvclassCount(eqvclassinputfile):
    txp_in_eqv = dict()
    with open(eqvclassinputfile) as f:
        txpcount = int(f.readline()[:-1])
        eqvcount = int(f.readline()[:-1])
        id_txp = dict()
        for i in range(txpcount):
            txp = f.readline()[:-1]
            id_txp[i] = txp
            txp_in_eqv[txp] = 0
        for i in range(eqvcount):
            line = f.readline()[:-1]
            ll = line.split()
            groupSize = int(ll[0])
            for j in range(1, groupSize+1):
                id = int(ll[j])
                txp_in_eqv[id_txp[id]] += 1
    with open(eqvclassinputfile + ".out", 'w') as f:
        for k in txp_in_eqv.keys():
            f.write(k + " " + str(txp_in_eqv[k]) + "\n")

if __name__ == "__main__":
    eqvclassinputfile = cliParser(sys.argv[1:])
    eqvclassCount(eqvclassinputfile)