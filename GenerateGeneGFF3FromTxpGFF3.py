import readline
import re

with open("./transcript_only.gff3") as f:
    rlt = [] # 18052 unique gene
    text = f.readline()[:-1]
    last_found = str()
    while text:
        try:
            found = re.search('FBgn(.+?);', text).group(1)
            ll = text.split()
            ll[2] = "exon"
            ll[8] = "Parent=FBgn" + found
            if last_found == found:
                rlt[-1][3] = str(min(int(rlt[-1][3]), int(ll[3])))
                rlt[-1][4] = str(max(int(rlt[-1][4]), int(ll[4])))
            else:
                rlt.append(ll)
            last_found = found
        except AttributeError:
            found = ''
        text = f.readline()[:-1]
    # print(rlt)
    with open("gene.gff3", 'w') as f_out:
        for i in range(len(rlt)):
            line = "\t".join(rlt[i]) + "\n"
            f_out.write(line)