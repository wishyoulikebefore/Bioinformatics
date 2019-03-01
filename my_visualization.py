# Created by zty on 2019/1/11
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
import os

def spearman_scatterplot(data,var_x,var_y,xmax=0.5,xmin=-0.5,ymax=0.5,ymin=-0.5):
    f, ax = plt.subplots(figsize=(10, 10))
    ax = sns.scatterplot(x=var_x, y=var_y, data=data)
    ax.hlines(y=0,xmin=xmin,xmax=xmax,colors="red",linestyles="dashed")
    ax.vlines(x=0, ymin=ymin, ymax=ymax, colors="red", linestyles="dashed")
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    return ax


def plot_sig(ax, xstart, xend, ystart, yend, sig):
    """
    思路：先画一横与两竖，在添加***

    参数都为列表形式，同时呈现不同数据之间比较的结果
    """
    for i in range(len(xstart)):
        ax.plot(np.full((2, 1), xstart[i]), np.arange(ystart[i], yend[i], yend[i] - ystart[i] - 0.1), color="k")
        ax.plot(np.full((2, 1), xend[i]), np.arange(ystart[i], yend[i], yend[i] - ystart[i] - 0.1), color="k")
        ax.plot(np.arange(xstart[i], xend[i] + 0.1, xend[i] - xstart[i]), np.full((2, 1), yend[i]), color="k")
        ax.annotate(r'%s' % (sig[i]), xy=((xend[i] + xstart[i]) / 2, yend[i] + 1), color="red")
    return ax


def plot_gene_scatter(rnaseq, wd, prefix, gene, xlabel=None, ylabel=None, sd=None):
    filter_ranseq = rnaseq.loc[[gene], :].T
    filter_ranseq = np.log2(filter_ranseq + 1)
    filter_ranseq.rename(columns={gene: "transcriptional level"}, inplace=True)
    sampleList = filter_ranseq.index
    record = []
    os.chdir(wd)
    for sample in sampleList:
        try:
            df = pd.read_table("%s/%s%s" % (sample, sample, prefix), index_col=0)
        except FileNotFoundError:
            record.append(np.nan)
            continue
        record.append(df.ix[gene]["Weighted methylation level"])
    methy_df = pd.DataFrame(np.array(record).reshape(-1, 1), index=sampleList, columns=["methylation level"])
    merge_df = filter_ranseq.join(methy_df).dropna(how="any")
    merge_df["sample"] = np.where(merge_df.index.str.contains("T"), "Tumor", "Normal")

    sns.set(style="white", font_scale=2)
    f, ax = plt.subplots(figsize=(10, 10))
    ax = sns.scatterplot(x="methylation level", y="transcriptional level", data=merge_df, hue="sample")
    if not xlabel:
        ax.set_xlabel(xlabel, fontdict={"size": 20, "weight": "bold"})
    if not ylabel:
        ax.set_ylabel(ylabel, fontdict={"size": 20, "weight": "bold"})
    ax.set_title(gene, fontdict={"size": 20, "weight": "bold"})
    if not sd:
        f.savefig("%s_scatterplot.png" % (gene), dpi=100, bbox_inches='tight')
    else:
        f.savefig("%s/%s_scatterplot.png" % (sd, gene), dpi=100, bbox_inches='tight')
