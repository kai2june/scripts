TLE, this contains bugs.
USE GenerateTxpFastaFromDecoy.sh INSTEAD.

# import sys, getopt
# import readline

# def cliParser(argv):
#     decoyfile = ''
#     referencetxpfile = ''
#     try:
#         opts, args = getopt.getopt(argv, "hd:r:", ["decoyfile", "referencetxpfile"])
#     except getopt.GetoptError:
#         print('python3 GenerateTxpFastaFromDecoy.py -d <decoyfile> -r <referencetxpfile>')
#         sys.exit(2)
#     for opt, arg in opts:
#         if opt == '-h':
#             print('python3 GenerateTxpFastaFromDecoy.py -d <decoyfile> -r <referencetxpfile>')
#             sys.exit()
#         elif opt in ("-d", "--decoyfile"):
#             decoyfile = arg
#         elif opt in ("-r", "--referencetxpfile"):
#             referencetxpfile = arg
#     print ("decoyfile ", decoyfile)
#     print ("referencetxpfile", referencetxpfile)
#     return decoyfile, referencetxpfile

# def SubsampleTxpFastaFromDecoy(decoyfile, referencetxpfile):
#     decoys = dict()
#     with open(decoyfile) as f:
#         line = f.readline()[:-1]
#         while line:
#             print(line)
#             decoys[line] = ""
#             line = f.readline()[:-1]
#     with open(referencetxpfile) as f:
#         cur_txpname = ""
#         line = f.readline()[:-1]
#         while line:
#             if line[0] == '>':
#                 if line[1:] in decoys:
#                     cur_txpname = line[1:]
#                 else:
#                     cur_txpname = ""
#                 continue
#             if cur_txpname:
#                 decoys[cur_txpname] = decoys[cur_txpname] + line
#             line = f.readline()[:-1]
#     with open(referencetxpfile+"_subsample.fasta", 'w') as f:
#         for k in decoys.keys():
#             f.write(">" + k + "\n")
#             f.write(decoys[k] + "\n")
#     return decoys

# if __name__ == "__main__":
#     decoyfile, referencetxpfile = cliParser(sys.argv[1:])
#     decoys = SubsampleTxpFastaFromDecoy(decoyfile, referencetxpfile)