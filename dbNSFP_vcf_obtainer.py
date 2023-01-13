import pandas as pd
import os


def vcf_obtainer(file):
    df1 = pd.read_csv("dbNSFP4.3c_variant.chr21.gz", usecols=[*range(0,4), *range(37,151)], sep="\t", dtype=str, compression='gzip', nrows=100000)
    
    for column in df1.columns:
        df1[column] = df1[column].str.replace(";", "&")
    
    df1["#CHROM"] = 'chr' + df1["#chr"].astype(str)
    df1["POS"] = df1["pos(1-based)"]
    df1["ID"] = str(".")
    df1["REF"] = df1["ref"]
    df1["ALT"] = df1["alt"]
    df1["QUAL"] = str(".")
    df1["FILTER"] = str(".")
    df1["INFO"]= "SIFT_score=" + df1["SIFT_score"] + ";SIFT_converted_rankscore=" + df1["SIFT_converted_rankscore"] + ";SIFT_pred=" + df1["SIFT_pred"] + ";SIFT4G_score=" + df1["SIFT4G_score"] + ";IFT4G_converted_rankscore=" + df1["SIFT4G_converted_rankscore"] + ";SIFT4G_pred=" + df1["SIFT4G_pred"] + ";LRT_converted_rankscore=" + df1["LRT_converted_rankscore"] + ";LRT_pred=" + df1["LRT_pred"] + ";LRT_Omega=" + df1["LRT_Omega"] + ";MutationTaster_score=" + df1["MutationTaster_score"] + ";MutationTaster_converted_rankscore=" + df1["MutationTaster_converted_rankscore"] + ";MutationTaster_pred=" + df1["MutationTaster_pred"] + ";MutationTaster_model=" + df1["MutationTaster_model"] + ";MutationTaster_AAE=" + df1["MutationTaster_AAE"] + ";MutationAssessor_score=" + df1["MutationAssessor_score"] + ";MutationAssessor_rankscore=" + df1["MutationAssessor_rankscore"] + ";FATHMM_score=" + df1["FATHMM_score"] + ";FATHMM_converted_rankscore=" + df1["FATHMM_converted_rankscore"] + ";FATHMM_pred=" + df1["FATHMM_pred"] + ";PROVEAN_score=" + df1["PROVEAN_score"] + ";PROVEAN_converted_rankscore=" + df1["PROVEAN_converted_rankscore"] + ";PROVEAN_pred=" + df1["PROVEAN_pred"] + ";MetaSVM_score=" + df1["MetaSVM_score"] + ";MetaSVM_rankscore=" + df1["MetaSVM_rankscore"] + ";MetaSVM_pred=" + df1["MetaSVM_pred"] + ";MetaLR_score=" + df1["MetaLR_score"] + ";MetaLR_rankscore=" + df1["MetaLR_rankscore"] + ";MetaLR_pred=" + df1["MetaLR_pred"] + ";Reliability_index=" + df1["Reliability_index"] + ";MetaRNN_score=" + df1["MetaRNN_score"] + ";MetaRNN_rankscore=" + df1["MetaRNN_rankscore"] + ";MetaRNN_pred=" + df1["MetaRNN_pred"] + \
    ";M-CAP_score=" + df1["M-CAP_score"] + \
    ";M-CAP_rankscore=" + df1["M-CAP_rankscore"] + \
    ";M-CAP_pred=" + df1["M-CAP_pred"] + \
    ";MutPred_score=" + df1["MutPred_score"] + \
    ";MutPred_rankscore=" + df1["MutPred_rankscore"] + \
    ";MutPred_protID=" + df1["MutPred_protID"] + \
    ";MutPred_AAchange=" + df1["MutPred_AAchange"] + \
    ";MutPred_Top5features=" + df1["MutPred_Top5features"] + \
    ";MVP_score=" + df1["MVP_score"] + \
    ";MVP_rankscore=" + df1["MVP_rankscore"] + \
    ";MPC_score=" + df1["MPC_score"] + \
    ";MPC_rankscore=" + df1["MPC_rankscore"] + \
    ";PrimateAI_score=" + df1["PrimateAI_score"] + \
    ";PrimateAI_rankscore=" + df1["PrimateAI_rankscore"] + \
    ";PrimateAI_pred=" + df1["PrimateAI_pred"] + \
    ";DEOGEN2_score=" + df1["DEOGEN2_score"] + \
    ";DEOGEN2_rankscore=" + df1["DEOGEN2_rankscore"] + \
    ";DEOGEN2_pred=" + df1["DEOGEN2_pred"] + \
    ";BayesDel_addAF_score=" + df1["BayesDel_addAF_score"] + \
    ";BayesDel_addAF_rankscore=" + df1["BayesDel_addAF_rankscore"] + \
    ";BayesDel_addAF_pred=" + df1["BayesDel_addAF_pred"] + \
    ";BayesDel_noAF_score=" + df1["BayesDel_noAF_score"] + \
    ";BayesDel_noAF_rankscore=" + df1["BayesDel_noAF_rankscore"] + \
    ";BayesDel_noAF_pred=" + df1["BayesDel_noAF_pred"] + \
    ";LIST-S2_score=" + df1["LIST-S2_score"] + \
    ";LIST-S2_rankscore=" + df1["LIST-S2_rankscore"] + \
    ";LIST-S2_pred=" + df1["LIST-S2_pred"] + \
    ";Aloft_Fraction_transcripts_affected=" + df1["Aloft_Fraction_transcripts_affected"] + \
    ";Aloft_prob_Tolerant=" + df1["Aloft_prob_Tolerant"] + \
    ";Aloft_prob_Recessive=" + df1["Aloft_prob_Recessive"] + \
    ";Aloft_prob_Dominant=" + df1["Aloft_prob_Dominant"] + \
    ";Aloft_pred=" + df1["Aloft_pred"] + \
    ";Aloft_Confidence=" + df1["Aloft_Confidence"] + \
    ";DANN_score=" + df1["DANN_score"] + \
    ";DANN_rankscore=" + df1["DANN_rankscore"] + \
    ";fathmm-MKL_coding_score=" + df1["fathmm-MKL_coding_score"] + \
    ";fathmm-MKL_coding_rankscore=" + df1["fathmm-MKL_coding_rankscore"] + \
    ";fathmm-MKL_coding_pred=" + df1["fathmm-MKL_coding_pred"] + \
    ";fathmm-MKL_coding_group=" + df1["fathmm-MKL_coding_group"] + \
    ";fathmm-XF_coding_score=" + df1["fathmm-XF_coding_score"] + \
    ";fathmm-XF_coding_rankscore=" + df1["fathmm-XF_coding_rankscore"] + \
    ";fathmm-XF_coding_pred=" + df1["fathmm-XF_coding_pred"] + \
    ";Eigen-raw_coding=" + df1["Eigen-raw_coding"] + \
    ";Eigen-raw_coding_rankscore=" + df1["Eigen-raw_coding_rankscore"] + \
    ";Eigen-PC-raw_coding=" + df1["Eigen-PC-raw_coding"] + \
    ";Eigen-PC-raw_coding_rankscore=" + df1["Eigen-PC-raw_coding_rankscore"] + \
    ";Eigen-PC-phred_coding=" + df1["Eigen-PC-phred_coding"] + \
    ";integrated_fitCons_score=" + df1["integrated_fitCons_score"] + \
    ";integrated_fitCons_rankscore=" + df1["integrated_fitCons_rankscore"] + \
    ";integrated_confidence_value=" + df1["integrated_confidence_value"] + \
    ";GM12878_fitCons_score=" + df1["GM12878_fitCons_score"] + \
    ";GM12878_fitCons_rankscore=" + df1["GM12878_fitCons_rankscore"] + \
    ";GM12878_confidence_value=" + df1["GM12878_confidence_value"] + \
    ";H1-hESC_fitCons_score=" + df1["H1-hESC_fitCons_score"] + \
    ";H1-hESC_fitCons_rankscore=" + df1["H1-hESC_fitCons_rankscore"] + \
    ";H1-hESC_confidence_value=" + df1["H1-hESC_confidence_value"] + \
    ";HUVEC_fitCons_score=" + df1["HUVEC_fitCons_score"] + \
    ";HUVEC_fitCons_rankscore=" + df1["HUVEC_fitCons_rankscore"] + \
    ";HUVEC_confidence_value=" + df1["HUVEC_confidence_value"] + \
    ";GERP++_NR=" + df1["GERP++_NR"] + \
    ";GERP++_RS=" + df1["GERP++_RS"] + \
    ";GERP++_RS_rankscore=" + df1["GERP++_RS_rankscore"] + \
    ";phyloP100way_vertebrate=" + df1["phyloP100way_vertebrate"] + \
    ";phyloP100way_vertebrate_rankscore=" + df1["phyloP100way_vertebrate_rankscore"] + \
    ";phyloP30way_mammalian=" + df1["phyloP30way_mammalian"] + \
    ";phyloP30way_mammalian_rankscore=" + df1["phyloP30way_mammalian_rankscore"] + \
    ";phyloP17way_primate=" + df1["phyloP17way_primate"] + \
    ";phyloP17way_primate_rankscore=" + df1["phyloP17way_primate_rankscore"] + \
    ";phastCons100way_vertebrate=" + df1["phastCons100way_vertebrate"] + \
    ";phastCons100way_vertebrate_rankscore=" + df1["phastCons100way_vertebrate_rankscore"] + \
    ";phastCons30way_mammalian=" + df1["phastCons30way_mammalian"] + \
    ";phastCons30way_mammalian_rankscore=" + df1["phastCons30way_mammalian_rankscore"] + \
    ";phastCons17way_primate=" + df1["phastCons17way_primate"] + \
    ";phastCons17way_primate_rankscore=" + df1["phastCons17way_primate_rankscore"] + \
    ";SiPhy_29way_pi=" + df1["SiPhy_29way_pi"] + \
    ";SiPhy_29way_logOdds=" + df1["SiPhy_29way_logOdds"] + \
    ";SiPhy_29way_logOdds_rankscore=" + df1["SiPhy_29way_logOdds_rankscore"] + \
    ";bStatistic=" + df1["bStatistic"] + \
    ";bStatistic_converted_rankscore=" + df1["bStatistic_converted_rankscore"] 
    df1.drop(df1.columns[:118], inplace=True, axis=1)
    df1['INFO'] = df1['INFO'].str.replace(" ", "_")
    df1 = df1.loc[df1["#CHROM"] != 'chr.' ]
    
        
    df2 = pd.read_csv("dbNSFP4.3c_variant.chr21.gz", usecols=[*range(2,4), *range(7,9), *range(37,151)], sep="\t", dtype=str, compression='gzip', nrows=100000)
    
    for column in df2.columns:
        df2[column] = df2[column].str.replace(";", "&")
        
    df2["#CHROM"] = 'chr' + df2["hg19_chr"].astype(str)
    df2["POS"] = df2["hg19_pos(1-based)"]
    df2["ID"] = str(".")
    df2["REF"] = df2["ref"]
    df2["ALT"] = df2["alt"]
    df2["QUAL"] = str(".")
    df2["FILTER"] = str(".")
    df2["INFO"]= "SIFT_score=" + df2["SIFT_score"] + ";SIFT_converted_rankscore=" + df2["SIFT_converted_rankscore"] + ";SIFT_pred=" + df2["SIFT_pred"] + ";SIFT4G_score=" + df2["SIFT4G_score"] + ";IFT4G_converted_rankscore=" + df2["SIFT4G_converted_rankscore"] + ";SIFT4G_pred=" + df2["SIFT4G_pred"] + ";LRT_converted_rankscore=" + df2["LRT_converted_rankscore"] + ";LRT_pred=" + df2["LRT_pred"] + ";LRT_Omega=" + df2["LRT_Omega"] + ";MutationTaster_score=" + df2["MutationTaster_score"] + ";MutationTaster_converted_rankscore=" + df2["MutationTaster_converted_rankscore"] + ";MutationTaster_pred=" + df2["MutationTaster_pred"] + ";MutationTaster_model=" + df2["MutationTaster_model"] + ";MutationTaster_AAE=" + df2["MutationTaster_AAE"] + ";MutationAssessor_score=" + df2["MutationAssessor_score"] + ";MutationAssessor_rankscore=" + df2["MutationAssessor_rankscore"] + ";FATHMM_score=" + df2["FATHMM_score"] + ";FATHMM_converted_rankscore=" + df2["FATHMM_converted_rankscore"] + ";FATHMM_pred=" + df2["FATHMM_pred"] + ";PROVEAN_score=" + df2["PROVEAN_score"] + ";PROVEAN_converted_rankscore=" + df2["PROVEAN_converted_rankscore"] + ";PROVEAN_pred=" + df2["PROVEAN_pred"] + ";MetaSVM_score=" + df2["MetaSVM_score"] + ";MetaSVM_rankscore=" + df2["MetaSVM_rankscore"] + ";MetaSVM_pred=" + df2["MetaSVM_pred"] + ";MetaLR_score=" + df2["MetaLR_score"] + ";MetaLR_rankscore=" + df2["MetaLR_rankscore"] + ";MetaLR_pred=" + df2["MetaLR_pred"] + ";Reliability_index=" + df2["Reliability_index"] + ";MetaRNN_score=" + df2["MetaRNN_score"] + ";MetaRNN_rankscore=" + df2["MetaRNN_rankscore"] + ";MetaRNN_pred=" + df2["MetaRNN_pred"] + \
    ";M-CAP_score=" + df2["M-CAP_score"] + \
    ";M-CAP_rankscore=" + df2["M-CAP_rankscore"] + \
    ";M-CAP_pred=" + df2["M-CAP_pred"] + \
    ";MutPred_score=" + df2["MutPred_score"] + \
    ";MutPred_rankscore=" + df2["MutPred_rankscore"] + \
    ";MutPred_protID=" + df2["MutPred_protID"] + \
    ";MutPred_AAchange=" + df2["MutPred_AAchange"] + \
    ";MutPred_Top5features=" + df2["MutPred_Top5features"] + \
    ";MVP_score=" + df2["MVP_score"] + \
    ";MVP_rankscore=" + df2["MVP_rankscore"] + \
    ";MPC_score=" + df2["MPC_score"] + \
    ";MPC_rankscore=" + df2["MPC_rankscore"] + \
    ";PrimateAI_score=" + df2["PrimateAI_score"] + \
    ";PrimateAI_rankscore=" + df2["PrimateAI_rankscore"] + \
    ";PrimateAI_pred=" + df2["PrimateAI_pred"] + \
    ";DEOGEN2_score=" + df2["DEOGEN2_score"] + \
    ";DEOGEN2_rankscore=" + df2["DEOGEN2_rankscore"] + \
    ";DEOGEN2_pred=" + df2["DEOGEN2_pred"] + \
    ";BayesDel_addAF_score=" + df2["BayesDel_addAF_score"] + \
    ";BayesDel_addAF_rankscore=" + df2["BayesDel_addAF_rankscore"] + \
    ";BayesDel_addAF_pred=" + df2["BayesDel_addAF_pred"] + \
    ";BayesDel_noAF_score=" + df2["BayesDel_noAF_score"] + \
    ";BayesDel_noAF_rankscore=" + df2["BayesDel_noAF_rankscore"] + \
    ";BayesDel_noAF_pred=" + df2["BayesDel_noAF_pred"] + \
    ";LIST-S2_score=" + df2["LIST-S2_score"] + \
    ";LIST-S2_rankscore=" + df2["LIST-S2_rankscore"] + \
    ";LIST-S2_pred=" + df2["LIST-S2_pred"] + \
    ";Aloft_Fraction_transcripts_affected=" + df2["Aloft_Fraction_transcripts_affected"] + \
    ";Aloft_prob_Tolerant=" + df2["Aloft_prob_Tolerant"] + \
    ";Aloft_prob_Recessive=" + df2["Aloft_prob_Recessive"] + \
    ";Aloft_prob_Dominant=" + df2["Aloft_prob_Dominant"] + \
    ";Aloft_pred=" + df2["Aloft_pred"] + \
    ";Aloft_Confidence=" + df2["Aloft_Confidence"] + \
    ";DANN_score=" + df2["DANN_score"] + \
    ";DANN_rankscore=" + df2["DANN_rankscore"] + \
    ";fathmm-MKL_coding_score=" + df2["fathmm-MKL_coding_score"] + \
    ";fathmm-MKL_coding_rankscore=" + df2["fathmm-MKL_coding_rankscore"] + \
    ";fathmm-MKL_coding_pred=" + df2["fathmm-MKL_coding_pred"] + \
    ";fathmm-MKL_coding_group=" + df2["fathmm-MKL_coding_group"] + \
    ";fathmm-XF_coding_score=" + df2["fathmm-XF_coding_score"] + \
    ";fathmm-XF_coding_rankscore=" + df2["fathmm-XF_coding_rankscore"] + \
    ";fathmm-XF_coding_pred=" + df2["fathmm-XF_coding_pred"] + \
    ";Eigen-raw_coding=" + df2["Eigen-raw_coding"] + \
    ";Eigen-raw_coding_rankscore=" + df2["Eigen-raw_coding_rankscore"] + \
    ";Eigen-PC-raw_coding=" + df2["Eigen-PC-raw_coding"] + \
    ";Eigen-PC-raw_coding_rankscore=" + df2["Eigen-PC-raw_coding_rankscore"] + \
    ";Eigen-PC-phred_coding=" + df2["Eigen-PC-phred_coding"] + \
    ";integrated_fitCons_score=" + df2["integrated_fitCons_score"] + \
    ";integrated_fitCons_rankscore=" + df2["integrated_fitCons_rankscore"] + \
    ";integrated_confidence_value=" + df2["integrated_confidence_value"] + \
    ";GM12878_fitCons_score=" + df2["GM12878_fitCons_score"] + \
    ";GM12878_fitCons_rankscore=" + df2["GM12878_fitCons_rankscore"] + \
    ";GM12878_confidence_value=" + df2["GM12878_confidence_value"] + \
    ";H1-hESC_fitCons_score=" + df2["H1-hESC_fitCons_score"] + \
    ";H1-hESC_fitCons_rankscore=" + df2["H1-hESC_fitCons_rankscore"] + \
    ";H1-hESC_confidence_value=" + df2["H1-hESC_confidence_value"] + \
    ";HUVEC_fitCons_score=" + df2["HUVEC_fitCons_score"] + \
    ";HUVEC_fitCons_rankscore=" + df2["HUVEC_fitCons_rankscore"] + \
    ";HUVEC_confidence_value=" + df2["HUVEC_confidence_value"] + \
    ";GERP++_NR=" + df2["GERP++_NR"] + \
    ";GERP++_RS=" + df2["GERP++_RS"] + \
    ";GERP++_RS_rankscore=" + df2["GERP++_RS_rankscore"] + \
    ";phyloP100way_vertebrate=" + df2["phyloP100way_vertebrate"] + \
    ";phyloP100way_vertebrate_rankscore=" + df2["phyloP100way_vertebrate_rankscore"] + \
    ";phyloP30way_mammalian=" + df2["phyloP30way_mammalian"] + \
    ";phyloP30way_mammalian_rankscore=" + df2["phyloP30way_mammalian_rankscore"] + \
    ";phyloP17way_primate=" + df2["phyloP17way_primate"] + \
    ";phyloP17way_primate_rankscore=" + df2["phyloP17way_primate_rankscore"] + \
    ";phastCons100way_vertebrate=" + df2["phastCons100way_vertebrate"] + \
    ";phastCons100way_vertebrate_rankscore=" + df2["phastCons100way_vertebrate_rankscore"] + \
    ";phastCons30way_mammalian=" + df2["phastCons30way_mammalian"] + \
    ";phastCons30way_mammalian_rankscore=" + df2["phastCons30way_mammalian_rankscore"] + \
    ";phastCons17way_primate=" + df2["phastCons17way_primate"] + \
    ";phastCons17way_primate_rankscore=" + df2["phastCons17way_primate_rankscore"] + \
    ";SiPhy_29way_pi=" + df2["SiPhy_29way_pi"] + \
    ";SiPhy_29way_logOdds=" + df2["SiPhy_29way_logOdds"] + \
    ";SiPhy_29way_logOdds_rankscore=" + df2["SiPhy_29way_logOdds_rankscore"] + \
    ";bStatistic=" + df2["bStatistic"] + \
    ";bStatistic_converted_rankscore=" + df2["bStatistic_converted_rankscore"] 
    
    df2.drop(df2.columns[:118], inplace=True, axis=1)
    df2['INFO'] = df2['INFO'].str.replace(" ", "_")
    df2 = df2.loc[df2["#CHROM"] != 'chr.' ]
    
    
    output_dir = os.getcwd()        
    output_file1 = os.path.join(output_dir, "hg38.vcf.gz" )
    df1.to_csv(output_file1, sep="\t", header=True, index=False)
    output_file2 = os.path.join(output_dir, "hg19.vcf.gz" )
    df2.to_csv(output_file2, sep="\t", header=True, index=False)
    


if __name__ == "__main__":
    files = [i for i in os.listdir(os.getcwd()) if i.endswith(".gz")]
    for file in files:
        vcf_obtainer(file) 
