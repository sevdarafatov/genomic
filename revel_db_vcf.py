import pandas as pd
import os

def vcf_obtainer(file):
    df1 = pd.read_csv(file, usecols=["chr" , "hg19_pos" , "ref" ,"alt", "aaref", "aaalt", "REVEL"] , sep=",", dtype=str)
    df1["#CHROM"] = 'chr' + df1["chr"].astype(str)
    df1["REF"] = df1["ref"]
    df1["ALT"] = df1["alt"]
    df1["POS"] = df1["hg19_pos"]
    df1["ID"] = str(".")
    df1["QUAL"] = str(".")
    df1["FILTER"] = str(".")
    df1["INFO"]= "aaref=" + df1["aaref"] + ";aaalt=" + df1["aaalt"] + ";REVEL=" + df1["REVEL"]
    df1= df1.drop(columns=["aaref", "aaalt", "REVEL"])
    df1 = df1[["#CHROM" , "POS" , "ID" ,"REF", "ALT", "QUAL", "FILTER", "INFO"]]



    df2 = pd.read_csv(file, usecols=["chr" , "grch38_pos" , "ref" ,"alt", "aaref", "aaalt", "REVEL"] , sep=",", dtype=str)
    df2["#CHROM"] = df2["chr"]
    df2["#CHROM"]='chr' + df2["chr"].astype(str)
    df2["REF"] = df2["ref"]
    df2["ALT"] = df2["alt"]
    df2["POS"] = df2["grch38_pos"]
    df2["ID"] = str(".")
    df2["QUAL"] = str(".")
    df2["FILTER"] = str(".")
    df2 = df2.loc[df2["grch38_pos"] != '.' ]
    df2["INFO"]= "aaref=" + df2["aaref"] + ";aaalt=" + df2["aaalt"] + ";REVEL=" + df2["REVEL"]
    df2 = df2.drop(columns=["aaref", "aaalt", "REVEL"])
    df2 = df2[["#CHROM", "POS" , "ID" ,"REF", "ALT", "QUAL", "FILTER", "INFO"]]


    output_dir = os.getcwd()        
    output_file1 = os.path.join(output_dir, "hg19.vcf" )
    output_file2 = os.path.join(output_dir, "GRCh38.vcf" )
    df1.to_csv(output_file1, sep="\t", header=True, index=False)
    df2.to_csv(output_file2, sep="\t", header=True, index=False)

if __name__ == "__main__":
    files = [i for i in os.listdir(os.getcwd()) if i.endswith(".csv")]
    for file in files:
        vcf_obtainer(file)
