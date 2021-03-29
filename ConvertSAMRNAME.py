import sys, getopt
import readline

def cliParser(argv):
    saminputfile = ''
    samoutputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hs:o:", ["saminputfile", "samoutputfile"])
    except getopt.GetoptError:
        print("python3 ConvertSAMRNAME.py -s <saminputfile> -o <samoutputfile>")
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("python3 ConvertSAMRNAME.py -s <saminputfile> -o <samoutputfile>")
            sys.exit()
        elif opt in ("-s", "--saminputfile"):
            saminputfile = arg
        elif opt in ("-o", "--samoutputfile"):
            samoutputfile = arg
    print("saminputfile ", saminputfile)
    print("samoutputfile ", samoutputfile)
    return saminputfile, samoutputfile

def convertRNAME(saminputfile, samoutputfile):
    out_vec = []
    with open(saminputfile) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            ll[2] = ll[0].split(":")[2]
            out_vec.append("\t".join(ll))
            line = f.readline()[:-1]
    with open(samoutputfile, 'w') as f:
        for i in out_vec:
            f.write(i + "\n")

if __name__ == "__main__":
    saminputfile, samoutputfile = cliParser(sys.argv[1:])
    convertRNAME(saminputfile, samoutputfile)