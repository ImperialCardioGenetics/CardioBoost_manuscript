import sys
from itertools import *
import re
from Bio import SeqIO
import pandas as pd
from collections import OrderedDict
from Bio import AlignIO
import Bio.SubsMat.MatrixInfo as SubsMat
import numpy as np

def find_transcript(disease):
    cantranscript={};
    if(disease=="arm"):
        gene_f="../../data/arrhythmia/arm_gene.txt"

    if(disease=="cm"):
        gene_f="../../data/cardiomyopathy/cm_gene.txt"

    gene=pd.read_table(gene_f,index_col=False)

    for i, row in gene.iterrows():
        symbol=row["gene_symbol"]
        transcript=row['gene_canon_transcript']
        cantranscript[symbol]=transcript

    return cantranscript;

#grantham matrix file
grantham_mat_file="/Users/xzhang13/Desktop/OneDrive/OneDrive - Imperial College London/Project1/LQTS/OtherScores/gratham_dict/grantham.matrix.txt"

#pam1 matrix file
pam1_mat_file="/Users/xzhang13/Desktop/OneDrive/OneDrive - Imperial College London/Project1/LQTS/OtherScores/pam1/pam1.txt"

aa_triple = {'Ala': 'A', 'Arg':'R', 'Asn':'N', 'Asp':'D', 'Cys':'C', 'Gln':'Q', 'Glu':'E', 'Gly':'G', 'His':'H', 'Ile':'I','Leu':'L', 'Lys':'K',
            'Met':'M', 'Phe':'F', 'Pro':'P', 'Pyl':'O', 'Ser':'S', 'Sec':'U', 'Thr':'T', 'Trp':'W', 'Tyr':'Y', 'Val':'V', 'Asx':'B', 'Glx':'Z',
            'Xaa':'X', 'Xle':'J'}

Gap_penality = -10

#subset of species

primate=["panTro4","ponAbe2","nomLeu3","rheMac3","papHam1","chlSab1","otoGar3","macFas5","gorGor3","calJac3","saiBol1"]
euar=["octDeg1","chiLan1","criGri1","tupChi1","mesAur1","cavPor3","jacJac1","mm10","hetGla2","ochPri3","micOch1","oryCun2","rn5","speTri2"]
laur=["vicPac2","camFer1","eptFus1","pteAle1","felCat5","bosTau7","myoDav1","canFam3","turTru2","capHir1","musFur1","eriEur2","equCab2","orcOrc1","pteVam1","myoLuc2","odoRosDiv1","ailMel1","susScr3","oviAri3","sorAra2","conCri1","panHod1","lepWed1","cerSim1"]
afro=["oryAfe1","eleEdw1","chrAsi1","loxAfr3","triMan1","echTel2"]
mam=["dasNov3","monDom5","ornAna1","sarHar1","macEug2"]
aves=["melUnd1","galGal4","ficAlb2","anaPla1","geoFor1","amaVit1","falPer1","colLiv1","falChe1","araMac1","pseHum1","melGal1","zonAlb1","taeGut2"]
sarc=["allMis1","pelSin1","latCha1","cheMyd1","anoCar2","chrPic1","apaSpi1","xenTro7"]
fish=["gadMor1","hapBur1","fr3","petMar2","oryLat2","astMex1","oreNil2","neoBri1","punNye1","xipMac1","lepOcu1","gasAcu1","tetNig2","takFla1","mayZeb1","danRer7"]



def Blosum62_score(ref_aa,alt_aa):
    value = SubsMat.blosum62.get((ref_aa,alt_aa))
    if value is None:
        value = SubsMat.blosum62.get((alt_aa,ref_aa))
    return(int(value))


def Pam250_score(ref_aa,alt_aa):
    value = SubsMat.pam250.get((ref_aa,alt_aa))
    if value is None:
        value = SubsMat.pam250.get((alt_aa,ref_aa))
    return(int(value))

###make_grantham_dict is copied from https://gist.github.com/arq5x/5408712
def make_grantham_dict(grantham_mat_file):
    f = open(grantham_mat_file)
    header = f.next().strip().split('\t')
    idx_to_aa = dict(zip(range(0,len(header)), header))

    grantham_dict = {}
    for line in f:
        fields = line.strip().split('\t')
        from_aa = fields[0]

        for idx, score in enumerate(fields):
            if idx == 0:
                continue
            to_aa = idx_to_aa[idx]
            grantham_dict[(from_aa, to_aa)] = score

    return grantham_dict



grantham_dict=make_grantham_dict(grantham_mat_file)


def Grantham_score(ref_aa,alt_aa):
    gscore=grantham_dict[(ref_aa,alt_aa)]
    return(int(gscore))


def make_pam1_dict(pam1_mat_file):
    f = open(pam1_mat_file)
    header = f.next().strip().split()
    idx_from_aa = dict(zip(range(0,len(header)), header))
    pam1_dict = {}
    for line in f:
        fields = line.strip().split()
        to_aa = fields[0]
        for idx, score in enumerate(fields):
            if idx == 0:
                continue
            from_aa = idx_from_aa[idx]
            pam1_dict[(from_aa, to_aa)] = score
            #print from_aa
    for idx,from_aa in idx_from_aa.items():
        if idx==0:
            continue
        pam1_dict[(from_aa,'Z')]=float(pam1_dict[(from_aa,'Q')])+float(pam1_dict[(from_aa,'E')])
        pam1_dict[(from_aa,'B')]=float(pam1_dict[(from_aa,'N')])+float(pam1_dict[(from_aa,'D')])
        pam1_dict[(from_aa,'X')]=Gap_penality
    return pam1_dict

pam1_dict=make_pam1_dict(pam1_mat_file)


def pam1_score(ref_aa,alt_aa):
    pscore=pam1_dict[(ref_aa,alt_aa)]
    return(int(pscore))


def extract_hgvsc(hgvsc):
    match=re.search('ENST\d+.\d+:c.(\d+)(\w)>(\w)',hgvsc)
    dict_nt={}
    dict_nt['pos']=match.group(1)
    dict_nt['ref']=match.group(2)
    dict_nt['alt']=match.group(3)
    return(dict_nt)

def extract_hgvsp(hgvsp):
    #ENSP00000433295.1:p.Arg41Cys
    match=re.search('ENSP\d+.\d:p.(\w\w\w)(\d+)(\w\w\w)',hgvsp)
    dict_aa={}
    dict_aa['ref']=match.group(1)
    dict_aa['pos']=match.group(2)
    dict_aa['alt']=match.group(3)
    return(dict_aa)


def msa_len(msa_file):
    record=SeqIO.parse(msa_file, "fasta").next()
    msa_length = len(record.seq)

    return(msa_length)

def blosum62_aa(msa_path,transcript,ref,alt,pos):
    pattern = re.compile(transcript)
    ref_dict=OrderedDict();
    alt_dict=OrderedDict();
    msa_file=open(msa_path,"rU");
    for record in SeqIO.parse(msa_file, "fasta"):
        if (pattern.match(record.id)):
            ref_id = record.id.replace(transcript, "ref_aa_blosum62")
            alt_id = record.id.replace(transcript, "alt_aa_blosum62")
            ref_dict[ref_id] = [int(Gap_penality if ref == "-" or record.seq[pos-1]=="-" else Blosum62_score(ref,record.seq[pos - 1]))]
            alt_dict[alt_id] = [int(Gap_penality if alt == "-" or record.seq[pos-1]=="-" else Blosum62_score(alt,record.seq[pos - 1]))]
    df1 = pd.DataFrame.from_dict(ref_dict).drop('ref_aa_blosum62_hg19', axis=1)
    #df1['ref_aa_sum'] = df1.sum(axis=1)
    df2 = pd.DataFrame.from_dict(alt_dict).drop('alt_aa_blosum62_hg19', axis=1)
    #df2['alt_aa_sum'] = df2.sum(axis=1)
    df = df1.join(df2)
    msa_file.close()
    #print list(df)
    return (df)


def pam1_aa(msa_path,transcript,ref,alt,pos):
    pattern = re.compile(transcript)
    ref_dict=OrderedDict();
    alt_dict=OrderedDict();
    msa_file=open(msa_path,"rU");
    for record in SeqIO.parse(msa_file, "fasta"):
        if (pattern.match(record.id)):
            ref_id = record.id.replace(transcript, "ref_aa_pam1")
            alt_id = record.id.replace(transcript, "alt_aa_pam1")
            ref_dict[ref_id] = [int(Gap_penality if ref == "-" or record.seq[pos-1]=="-" else pam1_score(ref,record.seq[pos - 1]))]
            alt_dict[alt_id] = [int(Gap_penality if alt == "-" or record.seq[pos-1]=="-" else pam1_score(alt,record.seq[pos - 1]))]
    df1 = pd.DataFrame.from_dict(ref_dict).drop('ref_aa_pam1_hg19', axis=1)
    #df1['ref_aa_sum'] = df1.sum(axis=1)
    df2 = pd.DataFrame.from_dict(alt_dict).drop('alt_aa_pam1_hg19', axis=1)
    #df2['alt_aa_sum'] = df2.sum(axis=1)
    df = df1.join(df2)
    msa_file.close()
    #print list(df)
    return (df)




def count_aa(msa_path,transcript,ref,alt,pos):
    pattern = re.compile(transcript)
    ref_dict=OrderedDict();
    alt_dict=OrderedDict();
    msa_file = open(msa_path, "rU");
    for record in SeqIO.parse(msa_file, "fasta"):
        if(pattern.match(record.id)):
            ref_id=record.id.replace(transcript, "ref_aa")
            alt_id=record.id.replace(transcript,"alt_aa")
#            ref_dict[ref_id]=[int(record.seq[pos-1]==ref)]
#            alt_dict[alt_id]=[int(record.seq[pos-1]==alt)]
            ref_dict[ref_id] = [np.nan if ref == "-" or record.seq[pos-1]=="-" else int(record.seq[pos-1]==ref)]
            alt_dict[alt_id] = [np.nan if alt == "-" or record.seq[pos-1]=="-" else int(record.seq[pos-1]==alt)]

    try:
        df1=pd.DataFrame.from_dict(ref_dict).drop('ref_aa_hg19', axis=1)
    except ValueError:
        print transcript
    #df1['ref_aa_sum']=df1.sum(axis=1)
    #print(df1['ref_aa_myoLuc2'])
    df2=pd.DataFrame.from_dict(alt_dict).drop('alt_aa_hg19', axis=1)
    #df2['alt_aa_sum']=df2.sum(axis=1)
    df_subset = count_aa_subset(df1, df2)
    msa_file.close()
    return(df_subset)


def ref_subset(subset_list):
    col_list=["ref_aa_" + subset_element for subset_element in subset_list]
    return(col_list)

def alt_subset(subset_list):
    col_list=["alt_aa_" + subset_element for subset_element in subset_list]
    return(col_list)

def count_aa_subset(df1,df2):
    df = pd.DataFrame()
    df['ref_aa_all']=df1.sum(axis=1)/df1.notnull().sum(axis=1)
    df['alt_aa_all']=df2.sum(axis=1)/df1.notnull().sum(axis=1)
    df['gap_aa_all']=df1.isnull().sum(axis=1)/len(df1.columns)
    df['ratio_ortho_all']=float(len(df1.columns))/99
    subset=[primate,euar,laur,afro,mam,aves,sarc,fish]
    subset_name = ["primate","euar","laur","afro","mam","aves","sarc","fish"]
    ref_col=list(df1.columns)
    alt_col=list(df2.columns)
    for i in range(0,len(subset)):
        ref_subset_exist=[x for x in ref_subset(subset[i]) if x in ref_col ]
        alt_subset_exist=[x for x in alt_subset(subset[i]) if x in alt_col ]
        df['ratio_ortho_'+subset_name[i]]=float(len(df1[ref_subset_exist].columns))/len(subset[i])
        #print df['ratio_ortho_'+subset_name[i]].values
        #print df1[ref_subset_exist].notnull().sum(axis=1);
        no_gap = df1[ref_subset_exist].notnull().sum(axis=1).values
        if len(df1[ref_subset_exist].columns)==0 or no_gap==0:
            df['gap_aa_' + subset_name[i]] = 1
            df['ref_aa_' + subset_name[i]] = 0
            df['alt_aa_' + subset_name[i]] = 0
            continue;
        df['ref_aa_'+subset_name[i]]=float(df1[ref_subset_exist].sum(axis=1))/df1[ref_subset_exist].notnull().sum(axis=1)
        df['alt_aa_'+subset_name[i]]=float(df2[alt_subset_exist].sum(axis=1))/df1[ref_subset_exist].notnull().sum(axis=1)
        df['gap_aa_'+subset_name[i]]=float(df1[ref_subset_exist].isnull().sum(axis=1))/len(df1[ref_subset_exist].columns)

    return(df)

#Average ref_aa matchness and gap_aa observedness over 99 speciess among 10 nearest AA
def count_aa_region(msa_path,transcript,pos,region_size=10):
    pattern = re.compile(transcript)
    ref_dict = OrderedDict();
    msa_file = open(msa_path, "rU");
    ref_aa_list=[];
    gap_aa_list=[];
    df=pd.DataFrame();
    for record in SeqIO.parse(msa_file, "fasta"):
        if (pattern.match(record.id)):
            ref_id = record.id.replace(transcript, "ref_aa")
            if ref_id == "ref_aa_hg19":
                ref_hg19=record.seq;
                ref_len=len(record.seq);
                continue;
                # insensitive to strand direction
                # pay attention to index 0/1
            if ref_len-pos<int(round(region_size*0.5)):# if pos is at the right end of the seq
                ortho_aa = count_average_region(ref_hg19,record.seq,2*pos-ref_len,ref_len)
            elif pos<=int(round(region_size*0.5)): # if pos is at the right beginning of the seq
                ortho_aa = count_average_region(ref_hg19,record.seq,1,2*pos-1)
            else:
                ortho_aa = count_average_region(ref_hg19,record.seq,pos-int(round(region_size*0.5)),pos+int(region_size-round(region_size*0.5)))
            ref_aa_list.extend(ortho_aa);
    ref_aa_region_mean = np.nan if len(ref_aa_list)==sum(np.isnan(ref_aa_list)) else np.mean([x for x in ref_aa_list if np.isnan(x)!=True])
    #ref_aa_region_std = np.nan if len(ref_aa_list)==sum(np.isnan(ref_aa_list)) else np.std([x for x in ref_aa_list if np.isnan(x)!=True])
    gap_aa_region_mean = np.mean(np.isnan(ref_aa_list))
    #gap_aa_region_std =
    msa_file.close()
    return ref_aa_region_mean,gap_aa_region_mean


def count_average_region(ref_seq,ortho_seq,start,end):
    ortho_aa=[];
    for i in range(start-1,end):
        ref_match=np.nan if ref_seq[i-1] == "-" or ortho_seq[i-1] == "-" else int(ortho_seq[i - 1] == ref_seq[i-1])
        ortho_aa.append(ref_match);
    return ortho_aa



def count_nt(msa_path,transcript,ref,alt,pos):
    pattern = re.compile(transcript)
    ref_dict=OrderedDict();
    alt_dict=OrderedDict();
    msa_file=open(msa_path,"rU")
    for record in SeqIO.parse(msa_file, "fasta"):
        if(pattern.match(record.id)):
            ref_id=record.id.replace(transcript, "ref_nt")
            alt_id=record.id.replace(transcript,"alt_nt")
            ref_dict[ref_id] = [np.nan if ref == "-" or record.seq[pos-1]=="-" else int(record.seq[pos-1]==ref)]
            alt_dict[alt_id] = [np.nan if alt == "-" or record.seq[pos-1]=="-" else int(record.seq[pos-1]==alt)]

    df1=pd.DataFrame.from_dict(ref_dict).drop('ref_nt_hg19', axis=1)
    df1['ref_nt_all']=df1.sum(axis=1)/df1.notnull().sum(axis=1)
    df2=pd.DataFrame.from_dict(alt_dict).drop('alt_nt_hg19', axis=1)
    df2['alt_nt_all']=df2.sum(axis=1)/df1.notnull().sum(axis=1)
    df=df1.join(df2)
    msa_file.close()
    return df[['ref_nt_all','alt_nt_all']]




def find_codon_pos(pos):
    #return a start pos and end pos of the codon having nt at Position pos
    #pos assumed not to be 0
    pos_codon=[];
    if pos%3==0:
        pos_codon=[pos-3,pos]
    elif pos%3==1:
        pos_codon=[pos-1,pos+2]
    elif pos%3==2:
        pos_codon=[pos-2,pos+1]
    return(pos_codon)

def count_codon(msa_path,transcript,ref,alt,pos):
    pattern = re.compile(transcript)
    ref_transcript = transcript+"_hg19"
    pattern_ref = re.compile(ref_transcript)
    ref_dict=OrderedDict();
    alt_dict=OrderedDict();
    codon_start=find_codon_pos(pos)[0]
    codon_end=find_codon_pos(pos)[1]
    msa_file=open(msa_path,"rU")
    for record in SeqIO.parse(msa_file,"fasta"):
        if(pattern.match(record.id)):
            ref_id=record.id.replace(transcript, "ref_codon")
            alt_id=record.id.replace(transcript,"alt_codon")
            if(pattern_ref.match(record.id)):
                ref_codon = str(record.seq[codon_start:codon_end])
                alt_seq = str(record.seq).replace(ref,alt,pos-1)
                alt_codon = alt_seq[codon_start:codon_end]
            ref_dict[ref_id]=[np.nan if ref_codon=="---" or record.seq[codon_start:codon_end]=="---" else int(str(record.seq[codon_start:codon_end])==ref_codon)]
            alt_dict[alt_id]=[np.nan if alt_codon=="---" or record.seq[codon_start:codon_end]=="---" else int(str(record.seq[codon_start:codon_end])==alt_codon)]
    df1=pd.DataFrame.from_dict(ref_dict).drop('ref_codon_hg19', axis=1)
    df1['ref_codon_all']=df1.sum(axis=1)/df1.notnull().sum(axis=1)
    df2=pd.DataFrame.from_dict(alt_dict).drop('alt_codon_hg19',axis=1)
    df2['alt_codon_all']=df2.sum(axis=1)/df1.notnull().sum(axis=1)
    df=df1.join(df2)
    msa_file.close()
    return df[['ref_codon_all','alt_codon_all']]

