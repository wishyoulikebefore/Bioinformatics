import json
import pandas
import argparse

"""
利用imgt数据中的hsa，统计 和10X的对比一下
#数据格式真奇怪，END=长度
V_gene
{
    "baseSequence": "file://human_V_IGH.fasta#IGHV1-2*02",
    "name": "IGHV1-2*02",
    "geneType": "V",
    "isFunctional": true,
    "chains": [ "IGH" ],
    "meta": {
      "comments": "X62106|IGHV1-2*02|Homo sapiens|F|V-REGION|163..458|296 nt|1| | | | |296+24=320| | |"
    },
    "anchorPoints": {
      "FR1Begin": 0,
      "CDR1Begin": 75,
      "FR2Begin": 99,
      "CDR2Begin": 150,
      "FR3Begin": 174,
      "CDR3Begin": 285,
      "VEnd": 296
    }
}
J_gene
{
    "baseSequence": "file://human_J_IGH.fasta#IGHJ2*01",
    "name": "IGHJ2*01",
    "geneType": "J",
    "isFunctional": true,
    "chains": [ "IGH" ],
    "meta": {
      "comments": "J00256|IGHJ2*01|Homo sapiens|F|J-REGION|932..984|53 nt|2| | | | |53+0=53| | |"
    },
    "anchorPoints": {
      "JBegin": 0,
      "FR4Begin": 22,    #CDR3 end
      "FR4End": 53
    }
}
D_gene
{
    "uri": "file://human_D_IGH.fasta#IGHD1-14*01",
    "range": {
      "from": 0,
      "to": 17
    },
    "sequence": "GGTATAACCGGAACCAC"
}
"""

def parse_args():
    parser = argparse.ArgumentParser(description="trim R2文件")
    parser.add_argument("-f",help="imgt数据")
    args = parser.parse_args()
    return args.d

def process_json(_file):
    data = json.load(open(_file))
    for item in data["genes"]:
        gene_type = item["geneType"]
        gene_name = item["name"]
        subtype = gene_name.split("*")[0]
        functional = item["isFunctional"]  #if the coding region has an open reading frame without stop codon
        chains = item["chains"]
        try:
            sequence = item["sequence"]
        except Exception as e:
            sequence = "None"
        if gene_type == "V":
            CDR3_Vgene_length = item["anchorPoints"]["VEnd"] - item["anchorPoints"]["CDR3Begin"]
            Vgene_length = item["anchorPoints"]["VEnd"]
        elif gene_type == "J":
            CDR3_Jgene_length = item["anchorPoints"]["FR4Begin"] - item["anchorPoints"]["JBegin"]
            Jgene_length = item["anchorPoints"]["FR4Begin"]
        elif gene_type == "D":
            Dgene_length = len(item["sequence"])
        else:
            pass




if __name__ == "__main__":
    imgt_file = parse_args()
    process_json(imgt_file)
