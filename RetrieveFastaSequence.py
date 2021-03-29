import sys, getopt

def cliParser(argv):
    txtinputfile = ''
    referenceinputfileone = ''
    referenceinputfiletwo = ''
    isHisat2 = False
    try:
        opts, args = getopt.getopt(argv, "hf:1:2:t", ["txtinputfile", "referenceinputfileone", "referenceinputfiletwo", "hisat2"])
    except getopt.GetoptError:
        print('python3 RetrieveFastaSequence.py -f <txtinputfile> -1 <referenceinputfileone> -2 <referenceinputfiletwo> [-t = --hisat2]')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 RetrieveFastaSequence.py -f <txtinputfile> -1 <referenceinputfileone> -2 <referenceinputfiletwo> [-t = --hisat2]')
            sys.exit()
        elif opt in ("-f", "--txtinputfile"):
            txtinputfile = arg
        elif opt in ("-1", "--referenceinputfileone"):
            referenceinputfileone = arg
        elif opt in ("-2", "--referenceinputfiletwo"):
            referenceinputfiletwo = arg
        elif opt in ("-t", "--hisat2"):
            isHisat2 = True
    print ("txtinputfile: ", txtinputfile)
    print ("referenceinputfileone: ", referenceinputfileone)
    print ("referenceinputfiletwo: ", referenceinputfiletwo)
    print ("isHisat2: ", isHisat2)
    return txtinputfile, referenceinputfileone, referenceinputfiletwo, isHisat2

def retrieveFastaSequence(txtinputfile, referenceinputfileone, referenceinputfiletwo):
    unmapped = dict()
    with open(txtinputfile) as f:
        line = f.readline()[:-1]
        while line:
            unmapped['>' + line.split()[0]] = ''
            line = f.readline()[:-1]
    with open(referenceinputfileone) as f:
        line = f.readline()[:-1]
        while line:
            seq = f.readline()[:-1]
            if line in unmapped.keys():
                unmapped[line] = seq
            line = f.readline()[:-1]
    with open(referenceinputfiletwo) as f:
        line = f.readline()[:-1]
        while line:
            seq = f.readline()[:-1]
            if line in unmapped.keys():
                unmapped[line] = seq
            line = f.readline()[:-1]
    with open("unmapped.fasta", 'w') as f:
        for k in unmapped.keys():
            f.write(k + "\n")
            f.write(unmapped[k] + "\n")

def filterFastaSequence(txtinputfile, referenceinputfileone, referenceinputfiletwo, isHisat2):
    unmapped = dict()
    with open(txtinputfile) as f:
        line = f.readline()[:-1]
        while line:
            if isHisat2:
                unmapped[line] = ''
                unmapped[line.split("/")[0] + "/2"] = ''
                line = f.readline()[:-1] # read redundant sequence
            else:
                unmapped['>' + line.split()[0]] = ''
                unmapped['>' + line.split()[0].split("/")[0] + "/2"] = ''
            line = f.readline()[:-1]

    llone = []
    with open(referenceinputfileone) as f:
        line = f.readline()[:-1]
        while line:
            seq = f.readline()[:-1]
            if not line in unmapped.keys():
                llone.append(line)
                llone.append(seq)
            line = f.readline()[:-1]
    lltwo = []
    with open(referenceinputfiletwo) as f:
        line = f.readline()[:-1]
        while line:
            seq = f.readline()[:-1]
            if not line in unmapped.keys():
                lltwo.append(line)
                lltwo.append(seq)
            line = f.readline()[:-1]

    with open("new_1.fasta", 'w') as f:
        for i in llone:
            f.write(i + "\n")
    with open("new_2.fasta", 'w') as f:
        for i in lltwo:
            f.write(i + "\n")

if __name__ == "__main__":
    txtinputfile, referenceinputfileone, referenceinputfiletwo, isHisat2 = cliParser(sys.argv[1:])
    # retrieveFastaSequence(txtinputfile, referenceinputfileone, referenceinputfiletwo)
    filterFastaSequence(txtinputfile, referenceinputfileone, referenceinputfiletwo, isHisat2)