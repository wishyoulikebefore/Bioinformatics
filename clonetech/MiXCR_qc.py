# Created by zty on 2019/1/17
import os
import pandas as pd
import argparse
import glob
import numpy as np

"""
先进行以下操作
ls *out | sed 's/.out//' | while read line ; do less ${line}.out | sed -n '/Report/,$p' | sed -n '2,20p' | sed 's/not TCR\/IG//' | sed 's/\(\?\)//' | sed 's/()//' | sed 's/(/:/' | sed 's/%)//' | sed 's/ :/:/g' | sed 's/: /:/g' | sed 's/,//' > ${line}.mixcr ; done
提取数据方法需要优化一下，利用sh脚本
"""

def parse_args():
    parser = argparse.ArgumentParser(description="提取fastqc信息")
    parser.add_argument("-wd",help="work dir")
    args = parser.parse_args()
    return args.wd

def process(wd):
    os.chdir(wd)
    count_data = []
    proportion_data = []
    sample_idList = []
    select_index = ["Successfully aligned reads", "Alignment failed no hits",
                    "Alignment failed because of absence of V hits",
                    "Alignment failed because of absence of J hits",
                    "No target with both V and J alignments",
                    "Overlapped and aligned", "TRA chains", "TRB chains"]
    for _file in glob.glob("*mixcr"):
        sample_id = _file.replace(".mixcr","")
        sample_idList.append(sample_id)
        data = pd.read_table(_file,sep=":",header=None,names=["item","count","proportion"],index_col=0)
        filtered_data = data.ix[select_index]
        count_data.append(filtered_data["count"].tolist())
        proportion_data.append(filtered_data["proportion"].tolist())
    count_matrix = pd.DataFrame(np.array(count_data).T,columns = sample_idList,index=select_index)
    proportion_matrix = pd.DataFrame(np.array(proportion_data).T, columns=sample_idList,index=select_index)
    count_matrix.to_csv("MiXCR_QC_count.csv")
    proportion_matrix.to_csv("MiXCR_QC_prorportion.csv")

if __name__ == "__main__":
    wd = parse_args()
    process(wd)
