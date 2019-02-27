import pandas as pd
import numpy as np
from scipy.stats import lognorm
import argparse
import os

def pre_process(file_dir):

    #将melt数据用pivot_table转换，替换index和columns名字
    matrix_df = pd.read_table(os.path.join(file_dir,"matrix.mtx"),sep=" ",skiprows=3,header=None,names=["gene","cell","count"])
    barcode_df = pd.read_table(os.path.join(file_dir,"barcodes.tsv"), header=None, names=["barcode"])
    barcodeDict = dict(zip(barcode_df.index+1,barcode_df["barcode"]))
    gene_df = pd.read_table(os.path.join(file_dir,"genes.tsv"),header=None, names=["ENSG","gene_name"])
    gene_df["ENSG_gene"] = gene_df["ENSG"] + "_" + gene_df["gene_name"]
    geneDict = dict(zip(gene_df.index+1,gene_df["ENSG_gene"]))
    matrix_df = matrix_df.pivot_table(["count"],index=["gene"],columns="cell")["count"]   #变成多层结构
    matrix_df.rename(index=geneDict,columns=barcodeDict,inplace=True)

    #根据min.cells = 3, min.genes = 200，percent.mito<0.05清洗数据
    """
    思路:
    min.cells 和 min.genes：分别对行和列计算非NA值数量（True+True=2,bool值的加法）
    """
    filtered_index = matrix_df.index[matrix_df.notnull().sum(axis=1)>=3]
    filtered_columns = matrix_df.columns[matrix_df.notnull().sum() >= 200]
    matrix_df = matrix_df.loc[filtered_index,filtered_columns]
    mito_df = matrix_df.ix[matrix_df.index.str.contains("MT-")]
    mito_percent_df = mito_df.sum()/matrix_df.sum()
    matrix_df = matrix_df[mito_percent_df[mito_percent_df < 0.05].index]
    return matrix_df

def log_normalization(matrix_df,scale_factor=10000):
    #log normalization   ln((value/colsum)*scale.factor)
    #计算每个样本各基因的相对表达量
    col_sum = matrix_df.sum()
    logNorm_df = np.log1p(matrix_df/col_sum*scale_factor)
    return logNorm_df    #后续使用

def find_variable_genes(logNorm_df,mean="ExpMean",dispersion="logVMR",x_low=0.1,x_high=8,y_low=1,y_high=np.inf,
                        bins=20,bin_method="equal_width",top_genes=1000,selection_method="mean_var"):
    """
    计算各基因表达水平的均值和离散程度
    此处有坑：NA值，Cpp中exp计算expm1;python中NA会略过，导致值不同，可以借鉴expm1，先fillna(0)，计算e值再-1
    mean:
    dispersion: the variance to mean ratio (VMR)
    计算dispersion的 Z-score
    """
    temp_df = np.exp(logNorm_df.fillna(0)) - 1
    expmean = temp_df.mean(axis=1)
    gene_mean = np.log1p(expmean)
    variance = np.power(temp_df.sub(expmean,axis=0), 2).sum(axis=1)/(len(logNorm_df.columns)-1)
    gene_dispersion = np.log(variance/expmean).replace(-np.inf,0)
    """
    如何将cut的分组结果从gene_mean移到gene_dispersion上
    注意这两个都是series(最终方法:先组合成一个df)
    """
    merge_df = pd.DataFrame([gene_mean, gene_dispersion], index=["gene_mean", "gene_dispersion"]).T
    if bin_method == "equal_width":
        merge_df["x_bin"] = pd.cut(merge_df["gene_mean"], bins, labels=range(bins))
    elif bin_method == "equal_frequency":
        merge_df["x_bin"] = pd.qcut(merge_df["gene_mean"], bins, labels=range(bins))
    else:
        raise SystemExit('bin_method输入无效')
    #用分厢法再计算一次
    merge_df2 = merge_df.pivot_table(index="x_bin", columns=merge_df.index, values="gene_dispersion")
    bin_mean = merge_df2.mean(axis=1)       #分组计算，需要忽略NA值
    bin_std = merge_df2.std(axis=1)
    merge_df["bin_mean"] = list(bin_mean.ix[merge_df["x_bin"]])     #index不同
    merge_df["bin_std"] = list(bin_std.ix[merge_df["x_bin"]])
    merge_df["scaled_dispersion"] = (merge_df["gene_dispersion"] - merge_df["bin_mean"]) / merge_df["bin_std"]
    mean_variance_df = merge_df[["gene_mean", "gene_dispersion","scaled_dispersion"]]
    pass_cutoff_df = mean_variance_df[(mean_variance_df["gene_mean"] > x_low) &
                                   (mean_variance_df["gene_mean"] < x_high) &
                                   (mean_variance_df["scaled_dispersion"] > y_low) &
                                   (mean_variance_df["scaled_dispersion"] < y_high)]
    if selection_method == "dispersion":
        sorted_mean_variance_df = mean_variance_df.sort_values(by="gene_dispersion", ascending=False, inplace=True)
        variable_genes = mean_variance_df.index[:top_genes]
    elif selection_method == "mean_var":
        variable_genes = pass_cutoff_df.index
    else:
        raise SystemExit('Invalid selection method')
    return variable_genes




matrix_df = pre_process("C:/Users/zuo_tianyu0101/Desktop/test")
logNorm_df = log_normalization(matrix_df)
