import sys, getopt
import readline

def cliParser(argv):
    decoyinputfile = ''
    fastainputfile = ''
    subfastaoutputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hd:f:o:", ["decoyinputfile", "fastainputfile", "subfastaoutputfile"])
    except getopt.GetoptError:
        print('python3 SubsampleTxpFastaFromDecoy.py -d <decoyinputfile> -f <fastainputfile> -o <subfastaoutputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 SubsampleTxpFastaFromDecoy.py -d <decoyinputfile> -f <fastainputfile> -o <subfastaoutputfile>')
            sys.exit()
        elif opt in ("-d", "--decoyinputfile"):
            decoyinputfile = arg
        elif opt in ("-f", "--fastainputfile"):
            fastainputfile = arg
        elif opt in ("-o", "--outputfile"):
            subfastaoutputfile = arg
    print ("decoyinputfile ", decoyinputfile)
    print ("fastainputfile ", fastainputfile)
    print ("subfastaoutputfile ", subfastaoutputfile)
    return decoyinputfile, fastainputfile, subfastaoutputfile

def subsample(decoyinputfile, fastainputfile, fastaoutputfile):
    txp = dict()
    with open(decoyinputfile) as f:
        line = f.readline()[:-1]
        while line:
            txp[line] = []
            line = f.readline()[:-1]

    found = False
    cur_txp = ''
    with open(fastainputfile) as f:
        line = f.readline()[:-1]
        while line:
            if line[0] == '>':
                found = False
                cur_txp = line.split()[0].split(".")[0].split(">")[1]
                if cur_txp in txp.keys():
                    found = True
            else:                        
                if found:
                    txp[cur_txp].append(line)
            line = f.readline()[:-1]
    
    with open(subfastaoutputfile, 'w') as f:
        for k in txp.keys():
            f.write(">" + k + "\n")
            for i in txp[k]:
                f.write(i + "\n")
    return txp

if __name__ == "__main__":
    decoyinputfile, fastainputfile, subfastaoutputfile = cliParser(sys.argv[1:])
    subsample(decoyinputfile, fastainputfile, subfastaoutputfile)