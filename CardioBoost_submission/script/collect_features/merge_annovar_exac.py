
# coding: utf-8

import pandas as pd
import sys

#usage: merge_annovar_exac.py annovar.txt exac.txt > out.txt
annovar_path = sys.argv[1]
exac_path = sys.argv[2]
out = sys.argv[3]

annovar=pd.read_table(annovar_path,index_col=False,low_memory=False)
annovar.head(n=5)

exac=pd.read_table(exac_path,low_memory=False)
exac.head(n=5)

#Check whether key is unique in each file
annovar['key']=annovar['Chr'].map(str).str.cat([annovar['Start'].map(str),annovar['Ref'],annovar['Alt']], sep='_')

def remove_duplicate_allele(df,key):
    seen = set()
    uniq = []
    duplicate = []
    n=0
    for x in df[key]:
        if x not in seen:
            uniq.append(x)
            seen.add(x)
        else:
            duplicate.append(n)
        n=n+1

    df = df.drop(df.index[duplicate])
    return df


annovar=remove_duplicate_allele(annovar,"key")


namesList = ['chrom','pos','ref','alt']
exac.columns=namesList+exac.columns[4:].tolist()
exac['key']=exac['chrom'].map(str).str.cat([exac['pos'].map(str),exac['ref'],exac['alt']], sep='_')

exac=remove_duplicate_allele(exac,"key")

result=pd.merge(annovar,exac,on='key')

result["AF_Adj"]=result["AC_Adj"]/result["AN_Adj"]

select_column=['chrom','pos','ref','alt','gene','pathogenic','AF_Adj','HGVSc','HGVSp','SIFT_score',
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
 'fathmm-MKL_coding_score',
 'MetaSVM_score',
 'MetaLR_score',
 'integrated_fitCons_score',
 'integrated_confidence_value',
 'GERP++_RS','phyloP100way_vertebrate',
 'phyloP20way_mammalian','phastCons100way_vertebrate',
 'phastCons20way_mammalian',
 'SiPhy_29way_logOdds','GenoCanyon_score',
 'Eigen_coding_or_noncoding', 'Eigen-raw','Eigen-PC-raw',
 'MCAP',
 'REVEL','PARASCOREzSCORE']

sim_df=result[select_column]

sim_df.columns=['CHROM','POS','REF','ALT','gene','pathogenic','AF_Adj','HGVSc','HGVSp','SIFT_score',
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
 'GERP','phyloP100way_vertebrate',
 'phyloP20way_mammalian','phastCons100way_vertebrate',
 'phastCons20way_mammalian',
 'SiPhy_29way_logOdds','GenoCanyon_score',
 'Eigen_coding_or_noncoding', 'Eigen_raw','Eigen_PC_raw',
 'MCAP',
 'REVEL','PARASCOREzSCORE']
sim_df.to_csv(out,sep="\t",index=False,header=True)
