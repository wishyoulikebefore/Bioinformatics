# Created by zty on 2018/10/26
import shutil
import os
import pandas as pd
import argparse

"""
遗留问题
Per base sequence content
warning if the difference between A and T, or G and C is greater than 10% in any position
fail if the difference between A and T, or G and C is greater than 20% in any position
"""

def parse_args():
    parser = argparse.ArgumentParser(description="提取fastqc信息")
    parser.add_argument("-wd",help="work dir")
    args = parser.parse_args()
    return args.wd

def process(wd):
    os.chdir(wd)
    data = []
    fileList = []
    columns = ["fail","warn","baseQ30","baseQ20","seqQ30","seqQ20"]
    for sample_id in glob.glob("*"):
        if os.path.isdir(sample_id):
            print("开始处理%s" %(sample_id))
            fileList.append("%s_R1" %(sample_id))
            fail_item, warning_item, baseQ30, baseQ20, seqQ30, seqQ20 = decompress_and_extract("%s/%s_combined_R1_fastqc.zip" %(sample_id,sample_id))
            data.append([fail_item,warning_item,baseQ30,baseQ20,seqQ30,seqQ20])
            fileList.append("%s_R2" % (sample_id))
            fail_item, warning_item, baseQ30, baseQ20, seqQ30, seqQ20 = decompress_and_extract("%s/%s_combined_R2_fastqc.zip" % (sample_id, sample_id))
            data.append([fail_item,warning_item,baseQ30,baseQ20,seqQ30,seqQ20])
            print("处理完成%s" %(sample_id))
        else:
            pass
    matrix = pd.DataFrame(data,columns=columns,index=fileList)
    matrix.to_csv("fastqc_outcome.csv")

def decompress_and_extract(compressed_file):
    file_name_piece = compressed_file.split("_")
    shutil.unpack_archive(compressed_file)
    decompressed_dir = compressed_file.split("/")[-1].replace(".zip", "")
    os.chdir(decompressed_dir)
    summary_df = pd.read_table("summary.txt",sep="\t",header=None,names=["result","indicator","fastq"],index_col=1)
    fail_item = summary_df[summary_df["result"]=="FAIL"].index.tolist()
    warning_item = summary_df[summary_df["result"]=="WARN"].index.tolist()
    infoDict = get_rowNum()
    baseQ30,baseQ20 = process_base_quality(infoDict["base_quality_start"], infoDict["base_quality_rowNum"])
    seqQ30,seqQ20 = process_sequence_quality(infoDict["sequence_quality_start"], infoDict["sequence_quality_rowNum"])
    os.chdir("../")
    shutil.rmtree(decompressed_dir,ignore_errors=True)
    return fail_item,warning_item,baseQ30,baseQ20,seqQ30,seqQ20

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
    Q20_df = base_df[base_df["Median"]>=20]
    Q20_location = sum(Q20_df.index.str.contains("-"))*4 + len(Q20_df.index)
    loc_Q20_proportion = Q20_location/read_length * 100
    return round(loc_Q30_proportion, 1),round(loc_Q20_proportion, 1)

def process_sequence_quality(start_row,row_num):
    file_name = "fastqc_data.txt"
    sequence_df = pd.read_table(file_name,skiprows=start_row,nrows=row_num,header=None,
                            names=["Quality","Count"])
    reads_num = sum(sequence_df["Count"])
    Q30_df = sequence_df[sequence_df["Quality"]>=30]
    Q30_reads = sum(Q30_df["Count"])
    seq_Q30_proportion = Q30_reads/reads_num * 100
    Q20_df = sequence_df[sequence_df["Quality"]>=20]
    Q20_reads = sum(Q20_df["Count"])
    seq_Q20_proportion = Q20_reads/reads_num * 100
    return round(seq_Q30_proportion, 1), round(seq_Q20_proportion, 1)

if  __name__ == "__main__":
    wd = parse_args()
    process(wd)
