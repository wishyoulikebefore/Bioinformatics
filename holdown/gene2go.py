import os
import pandas as pd
import numpy as np
import re
from collections import defaultdict

"""
用于RNAseq数据标化处理后的注释，用于后续的富集分析
"""

gtfDict = {}
with open("gencode.v19.annotation.gtf") as f:
    for line in f:
        if line.startswith("#"):
            continue 
        lineList = line.rstrip().split("\t")
        if lineList[2] == "gene":
            ENSG = lineList[8].split(";")[0].split(" ")[1].strip("\"").split(".")[0]
            for item in lineList[8].split(";"):
                 [k,v] = item.strip().split(" ")
                 if k=="gene_type":
                     gene_type = v.strip("\"")
                     break
        gtfDict[ENSG] = gene_type

gene2goDict = defaultdict(list)
with open("gene2go.20170325") as f:
    for line in f:
        if line.startswith("#"):
            continue
        lineList = line.rstrip().split("\t")
        gene_id = lineList[1]
        go_id = lineList[2]
        go_term = lineList[5]
        category = lineList[7]
        gene2goDict[gene_id].append(",".join([go_id,go_term,category]))
        
"""
希望输出
ENSG00000211973_IGHV1-69 Ensembl_gene_type go_annoList/NA
"""
output = open("record","a")
with open("Homo_sapiens.gene_info.20170325") as f:
    for line in f:
        if line.startswith("#"):
            continue
        lineList = line.rstrip().split("\t")
        gene_id = lineList[1]
        symbol = lineList[2]
        ensg = re.findall("Ensembl:ENS[A-Z]*G[0-9]+", lineList[5])
        for item in ensg:
            ENSG = item.split(":")[1]
            ENSG_symbol = "%s_%s" %(ENSG,symbol)
            gene2go = gene2goDict.get(gene_id)
            gene_type = gtfDict.get(ENSG)
            if gene2go:
                output.write("%s\t%s\t%s\t%s\n" %(ENSG_symbol,gene_id,gene_type,";".join(gene2go)))
            else:
                output.write("%s\t%s\t%s\t%s\n" %(ENSG_symbol,gene_id,gene_type,"None"))
