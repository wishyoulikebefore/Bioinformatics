import matplotlib
matplotlib.use('Agg')
import pandas as pd
import  matplotlib.pyplot as plt
import argparse
import numpy as np
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description="绘制RNAseq的火山图")
    parser.add_argument("-f")
    parser.add_argument("-n",type=int,required=False,default=10)  # 点亮上调或下调前多少的基因
    parser.add_argument("-p",required=False,default=0.05 )
    parser.add_argument("-fc",required=False,default=2)
    args = parser.parse_args()
    return args.f,args.p,args.fc,args.n

def add_text(ax,df,y_axis):
    if len(df) > 0:
        ax.text(1.02, y_axis, "Top %s Up-regulated Genes" %(len(df)), fontsize=15, transform=ax.transAxes)
        for n,(x,y,gene) in enumerate(zip(df["log2FC"],df["padj2"],df.index)):
            ax.text(x, y, n+1, ha='center', va='center', fontsize=8)
            ax.text(1.02, y_axis-0.04*(n+1), "%s  %s" %(n+1,gene), fontsize=10, transform=ax.transAxes)
    else:
        return

def draw_volcano(csv_file,padj,fc,num):
    expression_matrix = pd.read_csv(csv_file,index_col=0)
    f, ax = plt.subplots(figsize=(10, 8))
    expression_matrix["padj2"] = -np.log10(expression_matrix["padj"])
    expression_matrix.loc[(expression_matrix["sig"] == True) & (expression_matrix["log2FC"] > 0), "type"] = "up"
    expression_matrix.loc[(expression_matrix["sig"] == True) & (expression_matrix["log2FC"] < 0), "type"] = "down"
    expression_matrix.loc[(expression_matrix["sig"] == False), "type"] = "normal"
    expression_matrix["size"] = np.where(expression_matrix["sig"],8,4)
    up_df_ = expression_matrix[expression_matrix["type"] == "up"].sort_values(by="log2FC",ascending=False)[:num]
    down_df_ = expression_matrix[expression_matrix["type"] == "down"].sort_values(by="log2FC",ascending=True)[:num]
    xlim = max(max(up_df_["log2FC"].abs()),max(down_df_["log2FC"].abs()))*1.1
    ylim = max(expression_matrix["padj2"])*1.1
    cmap = {"normal":"lightgrey","up":"red","down":"mediumpurple"}     #擦，原来palette接受dict类型
    sns.scatterplot(x="log2FC", y="padj2", hue="type", size="size", palette=cmap, legend=False, data=expression_matrix)
    ax.set_xlim([-xlim,xlim])
    ax.set_ylim([0,ylim])
    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("none")
    ax.hlines(-np.log10(padj),-xlim,xlim,colors="grey",linestyles = "dashed")
    ax.vlines(np.log2(1/fc),0,ylim,colors="grey",linestyles = "dashed")
    ax.vlines(np.log2(fc),0,ylim,colors="grey",linestyles="dashed")
    ax.set_xlabel("log2FC")
    ax.set_ylabel("-log10(padj)")
    ax.set_title("volcano plot")
    add_text(ax,up_df_,1)
    add_text(ax,down_df_,0.5)
    f.savefig('volcano.pdf', dpi=100, bbox_inches='tight')

if __name__ == "__main__":
    csv_file, padj, fc, num = parse_args()
    draw_volcano(csv_file,padj,fc,num)
