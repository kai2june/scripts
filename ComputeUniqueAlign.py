import sys, getopt
import readline

def cliParser(argv):
    samtxtinputfile = ''
    toolinputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hs:t:", ["samtxtinputfile", "toolinputfile"])
    except getopt.GetoptError:
        print('python3 ComputeUniqueAlign.py -s <samtxtinputfile> -t <toolinputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 ComputeUniqueAlign.py -s <samtxtinputfile> -t <toolinputfile>')
            sys.exit()
        elif opt in ("-s", "--samtxtinputfile"):
            samtxtinputfile = arg
        elif opt in ("-t", "--toolinputfile"):
            toolinputfile = arg
    print ("samtxtinputfile ", samtxtinputfile)
    print ("toolinputfile ", toolinputfile)
    return samtxtinputfile, toolinputfile

if __name__ == "__main__":
    samtxtinputfile, toolinputfile = cliParser(sys.argv[1:])
    print(samtxtinputfile, toolinputfile)
    tooldd = dict()
    dd = dict()
    with open(toolinputfile) as f:
        line = f.readline()[:-1]
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            tooldd[ll[0]] = float(ll[4])
            dd[ll[0]] = 0.0
            line = f.readline()[:-1]

    with open(samtxtinputfile) as f:
        star = 0
        line = f.readline()[:-1]
        while line:
            if line in dd.keys():
                dd[line] += 1.0
            if line == "*":
                star += 1
            line = f.readline()[:-1]
    
    i = 0
    for k in tooldd.keys():
        if abs(dd[k]/2.0 - tooldd[k]) > 20.0:
            print("{}th {} groundcount:{} salmoncount:{}".format(i, k, dd[k]/2.0, tooldd[k]))
            i+=1
    print("align 0 time readcount:",star, "paircount:", star/2)