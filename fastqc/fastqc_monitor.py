# Created by zty on 2018/11/5
import shutil
import os
import pandas as pd
import argparse

"""
嵌入shell中，利用$?捕获异常信息
因此需要主动报错（raise error）
"""

def parse_args():
    parser = argparse.ArgumentParser(description="提取fastqc信息")
    parser.add_argument("-f",help="需要提取信息的文件")
    args = parser.parse_args()
    return args.f

def decompress_and_extract(compressed_file):
    shutil.unpack_archive(compressed_file)
    decompressed_dir = compressed_file.split("/")[-1].replace(".zip","")
    os.chdir(decompressed_dir)
    summary_df = pd.read_table("summary.txt",sep="\t",header=None,names=["result","indicator","fastq"],index_col=1)
    fail_item = summary_df[summary_df["result"]=="FAIL"].index.tolist()
    warning_item = summary_df[summary_df["result"]=="WARN"].index.tolist()
    print("Fail: %s" %(fail_item))
    print("Warn: %s" %(warning_item))
    infoDict = get_rowNum()
    warn = 0
    if "Per base sequence quality" in fail_item:
        process_base_quality(infoDict["base_quality_start"],infoDict["base_quality_rowNum"])
        warn = 1
    if "Per sequence quality scores" in fail_item:
        process_sequence_quality(infoDict["sequence_quality_start"], infoDict["sequence_quality_rowNum"])
        warn = 1
    os.chdir("../")
    shutil.rmtree(decompressed_dir,ignore_errors=True)
    if warn == 1:
        raise ValueError

def get_rowNum():
    file_name = "fastqc_data.txt"
    infoDict = {"base_quality_start":0,"base_quality_rowNum":0,"sequence_quality_start":0,"sequence_quality_rowNum":0}
    with open(file_name) as f:
        for nu,line in enumerate(f):
            if "Per base sequence quality" in line:
                infoDict["base_quality_start"] = nu + 2
            elif "Per sequence quality scores" in line:
                infoDict["sequence_quality_start"] = nu + 2
            elif "END_MODULE" in line:
                if infoDict["base_quality_start"] != 0 and infoDict["base_quality_rowNum"] == 0:
                    infoDict["base_quality_rowNum"] = nu - infoDict["base_quality_start"]
                elif infoDict["sequence_quality_start"] != 0 and infoDict["sequence_quality_rowNum"] == 0:
                    infoDict["sequence_quality_rowNum"] = nu - infoDict["sequence_quality_start"]
                else:
                    pass
    return infoDict

def process_base_quality(start_row,row_num):
    file_name = "fastqc_data.txt"
    base_df = pd.read_table(file_name,skiprows=start_row,nrows=row_num,header=None,
                            names=["Base","Mean","Median","Lower Quartile","Upper Quartile","10th Percentile","90th Percentile"],
                            index_col=0)
    read_length = int(base_df.index[-1])
    Q30_df = base_df[base_df["Median"]>=30]
    Q30_location = sum(Q30_df.index.str.contains("-"))*4 + len(Q30_df.index)
    loc_Q30_proportion = Q30_location/read_length * 100
    print("Proportion of loc which median quality >= 30: %s%%" % (round(loc_Q30_proportion, 1)))
    Q20_df = base_df[base_df["Median"]>=20]
    Q20_location = sum(Q20_df.index.str.contains("-"))*4 + len(Q20_df.index)
    loc_Q20_proportion = Q20_location/read_length * 100
    print("Proportion of loc which median quality >= 20: %s%%" % (round(loc_Q20_proportion, 1)))

def process_sequence_quality(start_row,row_num):
    file_name = "fastqc_data.txt"
    sequence_df = pd.read_table(file_name,skiprows=start_row,nrows=row_num,header=None,
                            names=["Quality","Count"])
    reads_num = sum(sequence_df["Count"])
    Q30_df = sequence_df[sequence_df["Quality"]>=30]
    Q30_reads = sum(Q30_df["Count"])
    seq_Q30_proportion = Q30_reads/reads_num * 100
    print("Proportion of seq which mean quality >= 30: %s%%" % (round(seq_Q30_proportion, 1)))
    Q20_df = sequence_df[sequence_df["Quality"]>=20]
    Q20_reads = sum(Q20_df["Count"])
    seq_Q20_proportion = Q20_reads/reads_num * 100
    print("Proportion of seq which mean quality >= 20: %s%%" % (round(seq_Q20_proportion, 1)))

if  __name__ == "__main__":
    file_name = parse_args()
    decompress_and_extract(file_name)
