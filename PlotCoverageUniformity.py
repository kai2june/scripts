import sys, getopt
import readline
import seaborn as sns
from matplotlib import pyplot as plt 
import numpy as np

def cliParser(argv):
    fastainputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hf:", ["fastainputfile"])
    except getopt.GetoptError:
        print('python3 PlotCoverageUniformity.py -f <fastainputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 PlotCoverageUniformity.py -f <fastainputfile>')
            sys.exit()
        elif opt in ("-f", "--fastainputfile"):
            fastainputfile = arg
    print("fastainputfile: ", fastainputfile)
    return fastainputfile

def parseFasta(fastainputfile):
    bin_count = 100
    txp_len = dict()
    txp_binsize= dict()
    txp_bin = dict()
    total_bin1K = [0] * (bin_count+1)
    total_bin1K5K = [0] * (bin_count+1)
    total_bin5K = [0] * (bin_count+1)
    total_bin = [0] * (bin_count+1)
    with open(fastainputfile) as f:
        line = f.readline()[:-1]
        while line:
            ll = line.split(":")
            if not ll[2] in txp_binsize.keys():
                txp_len[ll[2]] = int(ll[4])
                txp_binsize[ll[2]] = float(ll[4]) / (bin_count)
                txp_bin[ll[2]] = [0] * (bin_count+1)

            left_bin = float(ll[5]) / txp_binsize[ll[2]]
            right_bin = float(ll[6].split("/")[0]) / txp_binsize[ll[2]]
            if left_bin >= bin_count:
                left_bin = int(bin_count)
            if right_bin >= bin_count:
                right_bin = int(bin_count)
            for i in range(int(left_bin), int(right_bin)+1):
                txp_bin[ll[2]][i] += 1
                total_bin[i] += 1
                if txp_len[ll[2]] >=0 and txp_len[ll[2]] < 1000:
                    total_bin1K[i] += 1
                elif txp_len[ll[2]] >= 1000 and txp_len[ll[2]] < 5000:
                    total_bin1K5K[i] += 1
                elif txp_len[ll[2]] >= 5000:
                    total_bin5K[i] += 1
                else:
                    raise Exception("txp_len < 0, problematic.")

            line = f.readline()[:-1] # read useless sequence in this case
            line = f.readline()[:-1] # read readheader

    parseFasta_plot([total_bin, total_bin1K, total_bin1K5K, total_bin5K], bin_count)
    return txp_len
#     parseFasta_plot(total_bin, "0_inf_bp_txp")
#     parseFasta_plot(total_bin1K, "0_1000_bp_txp")
#     parseFasta_plot(total_bin1K5K, "1000_5000_bp_txp")
#     parseFasta_plot(total_bin5K, "5000_inf_bp_txp")
# def parseFasta_plot(total_bin, description):
#     plt.scatter(np.arange(len(total_bin)-1), total_bin[:-1])
#     # plt.scatter(len(total_bin), total_bin[-1])
#     plt.title(description + " coverage x=[0,99], x=100(polyA)")
#     plt.savefig(description + '_coverage_uniformity.png')

def parseFasta_plot(arr, bin_count):
    color = ['red', 'green', 'blue', 'black']
    cluster = []
    cluster_name = ['all_txp', 'txpLT1K', 'txp1K5K', 'txpGT5K']
    for i in range(len(arr)):
        cluster.append(plt.scatter(np.arange(bin_count), arr[i][:-1], color=color[i]))
    plt.legend(cluster, cluster_name, loc='upper left')
    plt.title("coverage_uniformity_10millionreads")
    plt.savefig("coverage")
    
if __name__ == "__main__":
    fastainputfile = cliParser(sys.argv[1:])
    txp_len = parseFasta(fastainputfile)