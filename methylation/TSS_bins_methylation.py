# Created by zty on 2019/1/30
"""
立flag式的比对：简单但很麻烦的切换
"""

import pandas as pd
import numpy as np
import os
import glob
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="TSS")
    parser.add_argument("-f",help="TSS.np2kb.txt")
    parser.add_argument("-rd",default="/ehpcdata/zuotianyu/bed",help="reference directory")
    parser.add_argument("-sd")
    args = parser.parse_args()
    return args.f,args.rd,args.sd

def process(meth,unmeth,threshold):
    if meth + unmeth < 4:
        return 0,0
    if meth < threshold:
        return 0,meth+unmeth
    else:
        return meth,meth+unmeth

def process_file(_file,bin_dict,chain_dict,blackList,sd):
    file_name = os.path.basename(_file)
    sample_id = file_name.replace(".TSS.np2kb.txt","")
    now_gene = ""
    bins = []
    now_bin = 0
    numerator = 0
    denominator = 0
    start = 0
    end = 0
    data = []

    with open(_file) as f:
        for nu, line in enumerate(f):
            lineList = line.rstrip().split()
            gene = lineList[9]
            loc = int(lineList[1])
            meth = int(lineList[3])
            unmeth = int(lineList[4])
            meth_type = lineList[5]
            threshold = int(lineList[7])
            chain = chain_dict.get(gene)
            if gene in blackList:
                continue
            if now_gene != gene:
                if now_gene == "":
                    pass
                else:
                    data.append([now_gene, "bin_%s" % (now_bin), np.divide(numerator, denominator)])
                now_gene = gene
                if chain == "+":
                    bins = list(range(1, 21))
                else:
                    bins = list(range(20, 0, -1))
                # 基于pop弹出
                now_bin = bins.pop(0)
                numerator = 0
                denominator = 0
                rangeList = bin_dict["%s_bin%s" % (now_gene, now_bin)]
                start = rangeList[0]
                end = rangeList[1]
            elif loc > end:
                data.append([now_gene, "bin_%s" % (now_bin), np.divide(numerator, denominator)])
                now_bin = bins.pop(0)
                rangeList = bin_dict["%s_bin%s" % (now_gene, now_bin)]
                start = rangeList[0]
                end = rangeList[1]
                numerator = 0
                denominator = 0
            else:
                pass
            if meth_type == "CG":
                add_numerator, add_denominator = process(meth, unmeth, threshold)
                numerator += add_numerator
                denominator += add_denominator
            else:
                pass
        data.append([now_gene, "bin_%s" % (now_bin), np.divide(numerator, denominator)])
        matrix = pd.DataFrame(data, columns=["gene", "bins", "methylation level"])
        pivot_matrix = matrix.pivot_table(["methylation level"], index=["gene"], columns=["bins"]).fillna(0)
        pivot_matrix["methylation level"].to_csv("%s/%s_tss_methylation_level.csv" %(sd,sample_id))

if __name__ == "__main__":
    _file,rd,sd = parse_args()

    bin_dict = defaultdict(list)
    chain_dict = {}
    with open("%s/200bp_knownCanonical_TSS_np2kb.bed2" %(rd)) as f:
        for line in f:
            lineList = line.rstrip().split("\t")
            bin_dict[lineList[5]].append(int(lineList[1]))
            bin_dict[lineList[5]].append(int(lineList[2]))
    f.close()

    with open("%s/knownCanonical_TSS_np2kb.bed2" %(rd)) as f:
        for line in f:
            lineList = line.rstrip().split("\t")
            chain_dict[lineList[5]] = lineList[4]
    f.close()

    blackList = [i.rstrip() for i in open("%s/blacklist" %(rd))]

    process_file(_file, bin_dict, chain_dict, blackList, sd)
