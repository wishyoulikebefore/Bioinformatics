# Created by zty on 2018/9/2
from functools import wraps
import time
import math
from itertools import combinations
import os
import pandas as pd
import numpy as np

chrList = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13',
           'chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']

def timer(func):
    @wraps(func)
    def wrapper(*args,**kwargs):
        start_time = time.time()
        func(*args,**kwargs)
        end_time = time.time()
        print("共耗时%s秒" % (round(end_time - start_time,2)))
    return wrapper

def split_bed(input_df, sd, sn, chrList=chrList):
    """
    sd:save directory
    sn:save name
    """
    if not os.path.exists(sd):
        os.mkdir(sd)
    if not os.path.exists(os.path.join(sd, sn)):
        os.mkdir(os.path.join(sd, sn))
        for _chr in chrList:
            try:
                data = input_df[input_df["chr"] == _chr]
                try:
                    data = data.sort_values(by="start")
                except:
                    data = data.sort_values(by="loc")
                data.to_csv("%s/%s/%s.csv" % (sd, sn, _chr))
            except KeyError:
                pass
    else:
        return

def split_data(input_df, wd, chrList=chrList):
    for _chr in chrList:
        try:
            data = input_df.loc[_chr, :]
            try:
                data = data.sort_values(by="start")
            except:
                data = data.sort_values(by="loc")
            data.to_csv("%s/temp_%s.csv" % (wd,_chr))
        except KeyError:
            pass

def cal_methylation_by_chr(_chr, data_sd, bed_sd, bed_sn, output):
    try:
        chr_data = pd.read_csv("%s/temp_%s.csv" % (data_sd, _chr), index_col=0)
        chr_bed = pd.read_csv("%s/%s/%s.csv" % (bed_sd, bed_sn, _chr), index_col=0)
    except FileNotFoundError:
        return
    for index, row in chr_bed.iterrows():
        ###此处有问题
        start = row["start"]
        end = row["end"]
        element, denominator, level = cal_methylation(chr_data, start, end)
        output.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (_chr,start,end,element,denominator,level))
    os.remove("%s/temp_%s.csv" % (data_sd, _chr))

def cal_methylation(df, start, end):
    filter_df = df[(df["loc"] >= start) & (df["loc"] <= end)]
    meth_df = filter_df[filter_df["pvalue"] < 0.05]
    try:
        ### 需要保留分子分母用于后续计算
        #level = meth_df["meth"].sum() / (filter_df["meth"].sum() + filter_df["unmeth"].sum())
        element = meth_df["meth"].sum()
        denominator = filter_df["meth"].sum() + filter_df["unmeth"].sum()
        level = element/denominator
        return element,denominator,level
    except ZeroDivisionError:
        return 0

def cal_region_methylation_level(bed_df,data_df,data_sd,data_sn,bed_sd,bed_sn,chrList=chrList):
    split_bed(bed_df,bed_sd,bed_sn,chrList)
    split_data(data_df, data_sd,chrList)
    output = open("%s/%s_DMR_region_methylation_level.txt" %(data_sd,data_sn), "a")
    output.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("chr", "start","end","element","denominator","Weighted methylation level"))
    for _chr in chrList:
        cal_methylation_by_chr(_chr, data_sd, bed_sd, bed_sn, output)

def normalize_rnaseq(rnaseq):
    select_df = rnaseq.ix[rnaseq.any(axis=1)]
    samples = rnaseq.columns
    f75 = rnaseq.apply(lambda x: np.percentile(x, 75))
    ref_column = np.argmin(abs(f75 - np.mean(f75)))
    normalized_factorDict = []
    for sample in samples:
        if sample == ref_column:
            normalized_factorDict.append(1)
        else:
            normalized_factorDict.append(cal_normalize_factor(select_df,sample,ref_column))
    sum_reads = select_df.sum()
    corrected_cpm = select_df * normalized_factorDict * 1000000 / sum_reads
    return corrected_cpm

def cal_normalize_factor(df, obs, ref, logratioTrim=0.3, sumTrim=0.05, cutoff=-1e10, weight=True):
    df = df[[obs,ref]]
    nO = df[obs].sum()
    nR = df[ref].sum()
    df = df.ix[df.all(axis=1)]
    # 牵扯到MA值计算中：-np.inf + n 的问题
    df["M"] = np.log2((df[obs] / nO) / (df[ref] / nR))
    df["A"] = (np.log2(df[obs] / nO) + np.log2(df[ref] / nR)) / 2
    df["var"] = (nO - df[obs]) / nO / df[obs] + (nR - df[ref]) / nR / df[ref]
    filter_df = df.ix[(df["M"] != np.inf) & (df["M"] != -np.inf) & (df["M"] != np.nan) &
                      (df["A"] != -np.inf) & (df["A"] > cutoff)]
    if filter_df["M"].abs().max() < 1e-6:
        print("%s数据有问题" %(obs))
        return
    # TMM normalization
    filter_genes_num = len(filter_df)
    low_L = np.floor(filter_genes_num * logratioTrim)
    high_L = filter_genes_num - low_L
    low_A = np.floor(filter_genes_num * sumTrim)
    high_A = filter_genes_num - low_A
    sort_M = filter_df["M"].sort_values(ascending=True).tolist()
    sort_A = filter_df["A"].sort_values(ascending=True).tolist()
    low_L_value, high_L_value = sort_M[int(low_L)], sort_M[int(high_L)]
    low_A_value, high_A_value = sort_A[int(low_A)], sort_A[int(high_A)]
    keep_df = filter_df.ix[(low_L_value <= filter_df["M"]) & (filter_df["M"] <= high_L_value) &
                           (low_A_value <= filter_df["A"]) & (filter_df["A"] <= high_A_value)]
    if weight:
        normalize_factor = np.sum(keep_df["M"] / keep_df["var"]) / np.sum(1 / keep_df["var"])
    else:
        normalize_factor = keep_df["M"].mean()
    return 2 ** normalize_factor

def slide_window(chr,start,end,chain,df,bin_length=500,bins=100):
    ### 默认针对076T.CX_report.txt.simplified
    ### print("本功能暂时只适用于CX_report.txt.simplified")
    if df.index.name == "chr":
        filter_df = df.ix[chr]
        filter_df = filter_df[(filter_df["loc"]>=start) & (filter_df["loc"]<=end)]
    else:
        filter_df = df[(df["chr"] == chr) & (df["loc"]>=start) & (df["loc"]<=end)]
    divisor, remainder = np.divmod(end-start+1-bin_length,bins-1)
    ###多出来的remainder加到前面的bin中
    sectionList = []
    if chain == "+":
        sectionList.append([start,start+bin_length-1])
        for i in range(1,remainder+1):
            each_add = divisor + 1
            sectionList.append([start + each_add * i, start + bin_length - 1 + each_add * i])
        for i in range(1,bins-remainder):
            each_add = divisor
            sectionList.append([start + each_add * i, start + bin_length - 1 + each_add * i])
    else:
        sectionList.append([end-bin_length+1, end])
        for i in range(1,remainder+1):
            each_add = divisor + 1
            sectionList.append([end - bin_length + 1 - each_add * i, end - each_add * i])
        for i in range(1,bins-remainder):
            each_add = divisor
            sectionList.append([end - bin_length + 1 - each_add * i, end - each_add * i])
    recordList = []
    for section in sectionList:
        a, b, level = cal_methylation(filter_df, section[0], section[1])
        recordList.append(level)
    return recordList
