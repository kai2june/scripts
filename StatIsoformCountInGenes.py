#python3 /home/scripts/StatIsoformCountInGenes.py -g /mammoth/flux_simulator_data/Drosophila_melanogaster.BDGP6.80.gff3 -p txp_1M.pro
import sys, getopt
import readline

def cliParser(argv):
    gff3inputfile = ''
    proinputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hg:p:", ["gff3inputfile", "proinputfile"])
    except getopt.GetoptError:
        print('python3 StatIsoformCountInGenes.py -g <gff3inputfile> -p <proinputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 StatIsoformCountInGenes.py -p <proinputfile>')
            sys.exit()
        elif opt in ("-g", "--gff3inputfile"):
            gff3inputfile = arg
        elif opt in ("-p", "--proinputfile"):
            proinputfile = arg
    print("gff3inputfile ", gff3inputfile)
    print("proinputfile ", proinputfile)
    return gff3inputfile, proinputfile

def statIsoformCountInGenes(gff3inputfile, proinputfile):
    gene_map_txp = dict()
    txp_map_count = dict()
    with open(gff3inputfile) as f:
        line = f.readline()[:-1]
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            if ll[2] == 'transcript':
                ll = ll[-1].split(";")
                gene = ll[1].split("=")[1]
                txp = ll[0].split("=")[1]
                if not gene in gene_map_txp.keys():
                    gene_map_txp[gene] = []
                gene_map_txp[gene].append(txp)
                txp_map_count[txp] = 0
            line = f.readline()[:-1]

    with open(proinputfile) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split()
            txp_map_count[ll[1]] = int(ll[9])
            line = f.readline()[:-1]

    gene_map_total = dict()
    gene_map_expressed = dict()
    for k in gene_map_txp.keys():
        if not len(gene_map_txp[k]) in gene_map_total:
            gene_map_total[len(gene_map_txp[k])] = 0
            gene_map_expressed[len(gene_map_txp[k])] = 0      
        gene_map_total[len(gene_map_txp[k])] += 1
        for k2 in gene_map_txp[k]:
            if txp_map_count[k2] > 10:
                gene_map_expressed[len(gene_map_txp[k])] += 1
                break
    for k in sorted(gene_map_total.keys()):
        print(k, gene_map_total[k], gene_map_expressed[k])
    
    decoy = []
    isoform_count_threshold = 1
    for k in gene_map_txp.keys():
        if len(gene_map_txp[k]) <= isoform_count_threshold:
            for txp in gene_map_txp[k]:
                decoy.append(txp)
    with open("txp_names_" + str(isoform_count_threshold) + ".txt", 'w') as f:
        for i in decoy:
            f.write(i + "\n")

if __name__ == "__main__":
    gff3inputfile, proinputfile = cliParser(sys.argv[1:])
    statIsoformCountInGenes(gff3inputfile, proinputfile)