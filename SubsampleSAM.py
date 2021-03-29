import sys, getopt
import readline
from enum import Enum
import numpy as np
import copy

### function
class SAM(Enum):
    QNAME = 0
    FLAG = 1
    RNAME = 2
    POS = 3
    MAPQ = 4
    CIGAR = 5
    RNEXT = 6
    PNEXT = 7
    TLEN = 8
    SEQ = 9
    QUAL = 10

def cliParser(argv):
    saminputfile = ''
    samoutputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hs:o:D:", ["saminputfile", "samoutputfile", "DEG"])
    except getopt.GetoptError:
        print('python3 SubsampleSAM.py -s <saminputfile> -o <samoutputfile> -D <DEG>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 SubsampleSAM.py -s <saminputfile> -o <samoutputfile> -D <DEG>')
            sys.exit()
        elif opt in ("-s", "--saminputfile"):
            saminputfile = arg
        elif opt in ("-o", "--samoutputfile"):
            samoutputfile = arg
        elif opt in ["-D", "--DEG"]:
            DEG = arg
    print("saminputfile", saminputfile)
    print("samoutputfile", samoutputfile)
    print("DEG", DEG)
    return saminputfile, samoutputfile, DEG

def parseSAM(saminputfile, samoutputfile, DEG):
    samrecords = []
    txpnames = set()
    with open(saminputfile) as f:
        line = f.readline()[:-1]
        if "SO:queryname" not in line:
            print("Exception: saminputfile should be sorted by queryname(i.e., samtools sort -n <saminputfile>)")
            sys.exit()
        while line[0] == '@':
            line = f.readline()[:-1]

        found = False
        cur_records = []
        cur_txpnames = set()
        while line:
            ll = line.split()
            if cur_records and (ll[SAM.QNAME.value] != cur_records[-1].split()[SAM.QNAME.value]):
                samrecords, txpnames, found, cur_records, cur_txpnames = parseSAM_found(samrecords, txpnames, found, cur_records, cur_txpnames)
            cur_records.append(line)
            cur_txpnames.add(ll[SAM.RNAME.value])

            if DEG in ll[SAM.RNAME.value]:
                found = True
            line = f.readline()[:-1]
        samrecords, txpnames, found, cur_records, cur_txpnames = parseSAM_found(samrecords, txpnames, found, cur_records, cur_txpnames)

    writeLines(samoutputfile + ".noheadersam", samrecords)
    writeLines(samoutputfile + ".txpnames", txpnames)
    return samrecords, txpnames

def parseSAM_found(samrecords, txpnames, found, cur_records, cur_txpnames):
    if found:
        for rec in cur_records:
            samrecords.append(rec)
        txpnames.update(cur_txpnames)
    cur_records.clear()
    cur_txpnames.clear()
    found = False
    return samrecords, txpnames, found, cur_records, cur_txpnames

def writeLines(filename, obj):
    with open(filename, 'w') as f_out:
        for i in obj:
            f_out.write(i + "\n")

if __name__ == "__main__":
    saminputfile, samoutputfile, DEG = cliParser(sys.argv[1:])
    samrecords, txpnames = parseSAM(saminputfile, samoutputfile, DEG)
