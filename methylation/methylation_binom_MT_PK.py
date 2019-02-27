# Created by zty on 2019/1/9
from scipy.stats import binom
import argparse
import subprocess
import pandas as pd
import sys
sys.path.append("/home/zuotianyu/script/stat")
from model_stats import binom_test, generate_binom_table
import multiprocessing
import glob
import os

"""
拆成24个文件，30min以内处理完

新增加功能：根据输入的conversion_rate先制作table

"""

def parse_args():
    parser = argparse.ArgumentParser(description="检测甲基化位点")
    parser.add_argument("-f")
    parser.add_argument("-fdr",required=False,default=0.01)   #non-conversion + sequencing error
    parser.add_argument("-sd", required=False, default="/ehpcdata/zuotianyu/WGBS_mid_file/binom_test")
    parser.add_argument("-pfx")     #输出文件前缀
    args = parser.parse_args()
    return args.f,args.fdr,args.sd,args.pfx

def get_EX(n,thresholdDict):
    exp = thresholdDict.get(n)
    if exp:
        return exp
    else:
        return 1

def get_pValue(mNum,n,fdr,binomDict):
    pValue = binomDict.get("%s_%s" % (n, mNum))
    if pValue:
        return pValue
    else:
        return binom_test(mNum,n,fdr)

def process_file(_file,fdr,binomDict,thresholdDict):

    output = open("%s.txt" %(_file),"a")
    with open(_file) as f:
        for line in f:
            lineList = line.rstrip().split()
            mNum = int(lineList[3])
            umNum = int(lineList[4])
            n = mNum+umNum
            EX = get_EX(n,thresholdDict)                    #默认0.05
            pValue = get_pValue(mNum,n,fdr,binomDict)
            output.write("%s\t%s\t%s\n" %("\t".join(lineList[:7]),EX,pValue))
    f.close()
    output.close()
    os.remove(_file)

if __name__ == "__main__":
    input_file,fdr,sd,output_prefix = parse_args()
    fdr = round(float(fdr),3)
    if not os.path.exists(output_prefix):
        os.mkdir(output_prefix)
    os.chdir(output_prefix)

    generate_binom_table(fdr,sd)
    table = pd.read_table("%s/binom_table_%s" %(sd,fdr), header=None, names=["con", "p"], dtype={"con": str, "p": str})
    table2 = pd.read_table("%s/binom_threshold_%s" %(sd,fdr), header=None, names=["n", "m"])
    binomDict = dict(zip(table["con"], table["p"]))
    thresholdDict = dict(zip(table2["n"], table2["m"]))

    wd = os.getcwd()
    output = open("%s.CX_report.txt" %(output_prefix), "a")
    subprocess.check_call(["split","-l","50000000",input_file,"%s_piece" %(output_prefix)])
    record = []
    pool = multiprocessing.Pool(processes=24)
    for piece_file in glob.glob("*piece*"):
        pool.apply_async(process_file, args=(piece_file,fdr,binomDict,thresholdDict,))
        record.append(piece_file)
    pool.close()
    pool.join()

    for item in record:
        _file = "%s.txt" %(item)
        with open(_file) as f:
            for line in f:
                output.write(line)
        f.close()
        os.remove(_file)
    output.close()
