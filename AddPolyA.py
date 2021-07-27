import readline
import sys, getopt

def cliParser(argv):
    txpfastainputfile = ''
    polyAlen = 0
    try:
        opts, args = getopt.getopt(argv, "ht:l:", ["help", "txpfastainputfile", "polyAlen"])
    except getopt.GetoptError:
        print('python3 AddPolyA.py -t <txpfastainputfile> -l <polyAlen>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 AddPolyA.py -t <txpfastainputfile>')
            sys.exit()
        elif opt in ("-t", "--txpfastainputfile"):
            txpfastainputfile = arg
        elif opt in ("-l", "--polyAlen"):
            polyAlen = int(arg)
    print ("txpfastainputfile ", txpfastainputfile)
    print ("polyAlen ", polyAlen)
    return txpfastainputfile, polyAlen

def parseTxpFasta(txpfastainputfile, polyAlen):
    txp_seq = dict()
    cur_txp = str()
    cur_seq = str()
    with open(txpfastainputfile) as f:
        line = f.readline()[:-1]
        while line:
            if line[0] == '>':
                if cur_seq:
                    cur_seq += 'A'*polyAlen
                    txp_seq[cur_txp] = cur_seq
                cur_txp = line
                cur_seq = str()
            else:
                cur_seq += line
            line = f.readline()[:-1]
        if cur_seq:
            cur_seq += 'A'*polyAlen
            txp_seq[cur_txp] = cur_seq
    for k in txp_seq.keys():
        print(k, txp_seq[k])

    outputFastaFileName = ".".join(txpfastainputfile.split('.')[:-1]) + ".polyA" + str(polyAlen) + ".fa"
    with open(outputFastaFileName, "w") as f:
        for k in txp_seq.keys():
            f.write(k + "\n")
            for i in range(len(txp_seq[k])):
                if i % 70 == 0 and i != 0:
                    f.write("\n")
                f.write(txp_seq[k][i])
            f.write("\n")

if __name__ == "__main__":
    txpfastainputfile, polyAlen = cliParser(sys.argv[1:])
    parseTxpFasta(txpfastainputfile, polyAlen)