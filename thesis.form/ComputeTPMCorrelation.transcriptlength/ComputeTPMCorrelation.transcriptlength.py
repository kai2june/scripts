import readline
import sys, getopt
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import operator

def cliParser(argv):
    transcript_length_inputfile = '' # quant.sf from Salmon
    transcript_tpm_inputfile = '' # TPM.txt from ComputeTPMCorrelation.py
    try:
        opts, args = getopt.getopt(argv, "hl:t:", ["transcript_length_inputfile", "transcript_tpm_inputfile"])
    except getopt.GetoptError:
        print('python3 ComputeTPMCorrelation.py -l <transcript_length_inputfile> -t <transcript_tpm_inputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 ComputeTPMCorrelation.py -l <transcript_length_inputfile> -t <transcript_tpm_inputfile>')
            sys.exit()
        elif opt in ("-l", "--transcript_length_inputfile"):
            transcript_length_inputfile = arg
        elif opt in ("-t", "--transcript_tpm_inputfile"):
            transcript_tpm_inputfile = arg

    print ("transcript_length_inputfile", transcript_length_inputfile)
    print ("transcript_tpm_inputfile", transcript_tpm_inputfile)
    return transcript_length_inputfile, transcript_tpm_inputfile

def getTranscriptLength(transcript_length_inputfile):
    transcript_length = dict()
    with open(transcript_length_inputfile, 'r') as f:
        line = f.readline()[:-1]
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            if ll[-1] == "NotInGroundFile" or ll[-1] == "NotInToolIndex":
                break
            transcript_length[ll[0]] = int(ll[1])
            line = f.readline()[:-1]
    return transcript_length

def getTranscriptTPM(transcript_tpm_inputfile):
    transcript_ground_tpm = dict()
    transcript_estimated_tpm = dict()
    with open(transcript_tpm_inputfile, 'r') as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            if ll[-1] == "NotInGroundFile" or ll[-1] == "NotInToolIndex":
                break
            transcript_ground_tpm[ll[0]] = float(ll[1])
            transcript_estimated_tpm[ll[0]] = float(ll[2])
            line = f.readline()[:-1]
    return transcript_ground_tpm, transcript_estimated_tpm

def computeTPMCorrelationGivenTranscriptLength(transcript_length, transcript_ground_tpm, transcript_estimated_tpm, lower_length, upper_length):
    ground = []
    estimated = []
    for k in transcript_estimated_tpm.keys():
        if (transcript_length[k] > lower_length and transcript_length[k] <= upper_length) or (upper_length == -1 and transcript_length[k] > lower_length):
            if not k in transcript_ground_tpm:
                exit(3)
            else:
                ground.append(transcript_ground_tpm[k])
            if not k in transcript_estimated_tpm:
                exit(4)
            else:
                estimated.append(transcript_estimated_tpm[k])
    spearman = str(stats.spearmanr(ground, estimated)).split("=")[1].split(",")[0]
    pearson = str(stats.pearsonr(ground, estimated)).split("(")[1].split(",")[0]
    # print("length = (%d, %d]" %(lower_length, upper_length))
    # print("transcript count: %d" %len(ground))
    filename = "corr_transcript_length_" + str(lower_length) + "_" + str(upper_length) + ".txt"
    with open(filename, 'w') as f:
        f.write("spearman {}\n".format(spearman))
        f.write("pearson {}\n".format(pearson))    

if __name__ == "__main__":
    transcript_length_inputfile, transcript_tpm_inputfile = cliParser(sys.argv[1:])
    transcript_length = getTranscriptLength(transcript_length_inputfile)
    transcript_ground_tpm, transcript_estimated_tpm = getTranscriptTPM(transcript_tpm_inputfile)
    length_interval = [0, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000, -1]
    computeTPMCorrelationGivenTranscriptLength(transcript_length, transcript_ground_tpm, transcript_estimated_tpm, 0, -1)
    for i in range(len(length_interval)-1):
        computeTPMCorrelationGivenTranscriptLength(transcript_length, transcript_ground_tpm, transcript_estimated_tpm, length_interval[i], length_interval[i+1])