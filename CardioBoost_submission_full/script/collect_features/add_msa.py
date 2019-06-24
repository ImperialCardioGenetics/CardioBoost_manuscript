'''
Created on 13 Apr 2017

@author: xzhang13
'''
import new_features as nf
import pandas as pd
import re
import sys
import argparse


p = argparse.ArgumentParser()
p.add_argument("-i", "--input_table", help="input.tsv", required=True)
p.add_argument("-o", "--output_table", help="output_tsv", required=True)
p.add_argument("-d","--disease",help="arm or cm",required=True)
p.add_argument("-m","--msa_path",help="path of msa alignment file",required=True)
args = p.parse_args()


dis=args.disease
file_path=args.input_table  #the input file, the file after annotation from Annovar and gnomAD
data_file=pd.read_table(file_path,index_col=False,low_memory=False)
output_path=args.output_table
msa_path=args.msa_path



output_df=pd.DataFrame()
msa_column=[];
type=dis;
for i,row in data_file.iterrows():
    df=pd.DataFrame()
    gene = row['gene']
    msa_aa=msa_path+type+"_gene_aa"
    msa_nt=msa_path+type+"_gene_nt"
    transcript = nf.find_transcript(dis)[gene]
    hgvsc = row['HGVSc']
    hgvsp = row['HGVSp']
    dict_nt= nf.extract_hgvsc(hgvsc)
    dict_aa= nf.extract_hgvsp(hgvsp)
    grantham_score= nf.Grantham_score(nf.aa_triple[dict_aa['ref']],nf.aa_triple[dict_aa['alt']])
    blosum62= nf.Blosum62_score(nf.aa_triple[dict_aa['ref']],nf.aa_triple[dict_aa['alt']])
    pam250=nf.Pam250_score(nf.aa_triple[dict_aa['ref']],nf.aa_triple[dict_aa['alt']])
    df_aa=nf.count_aa(msa_aa,transcript,nf.aa_triple[dict_aa['ref']],nf.aa_triple[dict_aa['alt']],int(dict_aa['pos']))
    df_nt=nf.count_nt(msa_nt,transcript,dict_nt['ref'],dict_nt['alt'],int(dict_nt['pos']))
    df_codon=nf.count_codon(msa_nt,transcript,dict_nt['ref'],dict_nt['alt'],int(dict_nt['pos']))
    ref_aa_region_mean,gap_aa_region_mean=nf.count_aa_region(msa_aa,transcript,int(dict_aa['pos']),region_size=10)
    df=df_aa.join(df_nt)
    df=df.join(df_codon)
    df['ref_aa_region_mean']=ref_aa_region_mean
    df['gap_aa_region_mean']=gap_aa_region_mean
    msa_df=df
    msa_column=msa_df.columns
    df['grantham_score']=grantham_score
    df['blosum62']=blosum62
    df['pam250']=pam250

    output_df=output_df.append(df)



msa_column=list(msa_column)



select_column = ['CHROM',
 'POS',
 'REF',
 'ALT',
 'gene',
# 'gene_prior','gene_domain','domain_EPV',
 'pathogenic',
  'AF_Adj','gnomAD_AF',
 'HGVSc',
 'HGVSp',
 'grantham_score','blosum62','pam250',
 'SIFT_score',
 'Polyphen2_HDIV_score',
 'Polyphen2_HVAR_score',
 'LRT_score',
 'MutationTaster_score',
 'MutationAssessor_score',
 'FATHMM_score',
 'PROVEAN_score',
 'VEST3_score',
 'CADD_raw',
 'DANN_score',
 'fathmm_MKL_coding_score',
 'MetaSVM_score',
 'MetaLR_score',
 'integrated_fitCons_score',
 'integrated_confidence_value',
 'GERP', 'phyloP100way_vertebrate',
 'phyloP20way_mammalian', 'phastCons100way_vertebrate',
 'phastCons20way_mammalian',
 'SiPhy_29way_logOdds', 'GenoCanyon_score',
 'Eigen_coding_or_noncoding', 'Eigen_raw', 'Eigen_PC_raw',
 'MCAP',
 'REVEL','PARASCOREzSCORE','mis_badness','MPC']+msa_column

output_df=output_df.reset_index(drop=True)
output=data_file.join(output_df)
output=output[select_column]
output.to_csv(output_path,sep="\t",index=False,header=True)
