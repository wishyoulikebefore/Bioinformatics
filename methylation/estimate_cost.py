import pandas as pd
from dateparser import parse

"""
原始数据
ls -lh | awk '$5~"G"{print $9 "\t" $5 "\t" $6 " " $7 " " $8}' | sed '1d'
CHG_OT_R18046368LF01-BDESCCMET217-N_R1_val_1_bismark_bt2_pe.txt 107G    Oct 12 12:04
CHH_OB_R18046368LF01-BDESCCMET217-N_R1_val_1_bismark_bt2_pe.txt 392G    Oct 12 12:04
CHH_OT_R18046368LF01-BDESCCMET217-N_R1_val_1_bismark_bt2_pe.txt 397G    Oct 12 12:04
CpG_OB_R18046368LF01-BDESCCMET217-N_R1_val_1_bismark_bt2_pe.txt 23G     Oct 12 12:04
CpG_OT_R18046368LF01-BDESCCMET217-N_R1_val_1_bismark_bt2_pe.txt 23G     Oct 12 12:04
R18046368LF01-BDESCCMET217-N_R1_val_1_bismark_bt2_pe.bam        97G     Oct 12 00:00
R18046368LF01-BDESCCMET217-N_R1_val_1_bismark_bt2_pe.bedGraph.gz        6.0G    Oct 15 02:08
R18046368LF01-BDESCCMET217-N_R1_val_1_bismark_bt2_pe.bismark.cov.gz     5.7G    Oct 15 02:08
R18046368LF01-BDESCCMET217-N_R1_val_1_bismark_bt2_pe.CX_report.txt      33G     Oct 15 06:52
R18046368LF01-BDESCCMET217-N_R1_val_1.fq.gz     25G     Oct 10 02:53
R18046368LF01-BDESCCMET217-N_R2_val_2.fq.gz     25G     Oct 10 02:53
开始时间
Tue Oct 9 13:25:39 CST 2018 trim_galore.bismark.sh started...
后续计划：把Tumor数据加进来，再计算一下
"""

def rough_estimate(record_file,start_time):
    df = pd.read_table(record_file,header=None,names=["file_name","size","date"],index_col=0)
    df["standard_date"] = df.apply(lambda x:parse(x["date"]),axis=1)     #利用dateparser.parse将日期数据处理成标准写法
    df["standard_date"] = df.apply(lambda x:pd.Period(x["standard_date"],freq="H"),axis=1)
    df["start"] = pd.Period(parse(start_time),freq="H")
    df.drop(["date"],axis=1,inplace=True)
    df["generate_time"] = df["standard_date"] - df["start"]
    df["exist_time"] = max(df["generate_time"]) - df["generate_time"]
    df["size"] = df["size"].str.replace("G","")
    df["save_cost"] = df.apply(lambda x:float(x["size"])*x["exist_time"]*2.6/1000,axis=1)
    intermediate_save_cost = df["save_cost"].sum()
    fastq_save_cost = 69/1000*2.6*max(df["generate_time"])
    sum_save_cost = intermediate_save_cost + fastq_save_cost
    return sum_save_cost

def rough_estimate2(record_file,start_time):
    """
    改进之处：R18046368LF01-BDESCCMET217-N_R1_val_1_bismark_bt2_pe.CX_report.txt会以中间文件的形式存在
    以Tumor为例，中间文件最终汇总是1020G（约1个T）
    考虑N的测序量为T一半，因此预估0.5T，是最终文件的15倍
    不太了解具体的时间分布，开始是小文件，后面所有文件会一直存在一段时间，所以预估10倍
    至于这个文件的开始时间，应该是在val_1_bismark_bt2_pe.txt这类文件生成之后
    """
    df = pd.read_table(record_file,header=None,names=["file_name","size","date"],index_col=0)
    df["standard_date"] = df.apply(lambda x:parse(x["date"]),axis=1)
    df["standard_date"] = df.apply(lambda x:pd.Period(x["standard_date"],freq="H"),axis=1)
    df["start"] = pd.Period(parse(start_time),freq="H")
    df.drop(["date"],axis=1,inplace=True)
    df["generate_time"] = df["standard_date"] - df["start"]
    df["exist_time"] = max(df["generate_time"]) - df["generate_time"]
    df["size"] = df["size"].str.replace("G","")
    df["save_cost"] = df.apply(lambda x:float(x["size"])*x["exist_time"]*2.6/1000,axis=1)
    intermediate_save_cost = df["save_cost"].sum()
    fastq_save_cost = 69/1000*2.6*max(df["generate_time"])
    temp_df = df[df.index.str.contains("pe.CX_report.txt")]
    temp_start_time = max(df[df.index.str.contains("bismark_bt2_pe.txt")]["generate_time"])
    temp_cost = float(temp_df["size"])*(float(temp_df["generate_time"])-temp_start_time)*10*2.6/1000
    sum_save_cost = intermediate_save_cost + fastq_save_cost + temp_cost
    return sum_save_cost

rough_estimate2("file_data","Oct 9 13:25:39")
