import sys, getopt
import readline

def cliParser(argv):
    proinputfile = ''
    saminputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hg:s:", ["proinputfile", "saminputfile"])
    except getopt.GetoptError:
        print('python3 SAMFalsePositive.py <-g proinputfile> <-s saminputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 SAMFalsePositive.py <-g proinputfile> <-s saminputfile>')
            sys.exit()
        elif opt in ("-g", "--proinputfile"):
            proinputfile = arg
        elif opt in ("-s", "--saminputfileone"):
            saminputfile = arg
    print ("proinputfile: ", proinputfile)
    print ("saminputfile: ", saminputfile)
    return proinputfile, saminputfile

def countProCount(proinputfile):
    proCount = dict()
    with open(proinputfile) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            proCount[ll[1]] = int(ll[9])
            line = f.readline()[:-1]
    return proCount

def countSAMCount(saminputfile):
    samCount = dict()
    with open(saminputfile) as f:
        line = f.readline()[:-1]
        while line:
            if line[0] == '@':
                line = f.readline()[:-1]
                continue
            else:
                txp_name = line.split()[2]
                if not txp_name in samCount.keys():
                    samCount[txp_name] = 1
                else:
                    samCount[txp_name] += 1
            line = f.readline()[:-1]
    return samCount

if __name__ == "__main__":
    proinputfile, saminputfile = cliParser(sys.argv[1:])
    proCount = countProCount(proinputfile)
    samCount = countSAMCount(saminputfile)
    
    FP = []
    for k in samCount.keys():
        print("txp_name, samCount, proCount:", k, samCount[k], proCount[k])
        if samCount[k] > 0 and proCount[k] == 0:
            FP.append(k)
    with open("FP.txt", 'w') as f:
        for k in FP: 
            f.write("txp_name {}, samCount {}, proCount {}\n".format(k, samCount[k], proCount[k]))