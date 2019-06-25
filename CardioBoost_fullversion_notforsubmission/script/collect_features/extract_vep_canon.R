#!/usr/bin/env Rscript
library("optparse")
rm(list=ls())
require(dplyr)

option_list = list(
   make_option(c("-p", "--pathogenic"), type="character", default=NULL,
               help="with clinical significance", metavar="character"),
  make_option(c("-t","--transcript"),type="character",default=NULL,help="gene-canonical transcript file",metavar="character"),
  make_option(c("-v","--vep"),type="character",default=NULL,help="variants list with vep annotations",metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



#input file 1: variant list with clinvar interpretation
p_file=opt$pathogenic

#input file 2: vep_annotated file
vep_file=opt$vep;
#input file 2: canon transcript info file
##TODO: fixed format of gene infor file
gene_info_file=opt$transcript

#output file
output_file=opt$out


patho<-read.table(p_file,header=TRUE,sep="\t")
colnames(patho)<-tolower(colnames(patho))
patho[,"variation"]=paste(patho[,"chrom"],patho[,"pos"],patho[,"ref"],patho[,"alt"],sep = "_")
gene_info<-read.table(gene_info_file,header=TRUE,sep="\t")

#clinvar$variation<-paste(paste(clinvar[,'chrom'],clinvar[,'pos'],sep="_"),paste(clinvar[,'ref'],clinvar[,'alt'],sep="/"),sep="_")
vep_variants<-read.table(vep_file,header=FALSE,comment.char = "#",stringsAsFactors = FALSE)
colnames(vep_variants)<-c("Uploaded_variation","Location",	"Allele","Gene","Feature","Feature_type",
                           "Consequence",	"cDNA_position",	"CDS_position",	"Protein_position",	"Amino_acids",
                           "Codons","Existing_variation","ALLELE_NUM","IMPACT","DISTANCE",	"STRAND",	"FLAGS",	"MINIMISED",	"SYMBOL",
                            "SYMBOL_SOURCE",	"HGNC_ID",	"HGVSc",	"HGVSp",	"HGVS_OFFSET")

#select the columns or creat the columns needed: chr pos ref alt hgvsc hgvsp gene_symbol molecular_consequence(double check)







vep_can<-NULL
for (i in 1:nrow(gene_info)){
  gene_info_row<-gene_info[i,]
  gene_v<-subset(vep_variants,Feature==as.character(gene_info_row$gene_canon_transcript))
  vep_can<-rbind(vep_can,gene_v)
}

l=strsplit(as.character(vep_can$Uploaded_variation),"_")
df=data.frame(t(sapply(l,c)))
chr=df[,1]
pos=df[,2]
ref=df[,3]
alt=df[,4]

vep_selected=data.frame(variation=vep_can$Uploaded_variation,chrom=chr,pos=pos,ref=ref,alt=alt,gene=vep_can$SYMBOL,HGVSc=vep_can$HGVSc,HGVSp=vep_can$HGVSp)
#already check whether all variants are missense variants


df=merge(vep_selected,patho[,c("variation","pathogenic")],by="variation")[,-1]



frameshift=as.integer(grepl('Met1\\?',df$HGVSp))
df=df[frameshift!=1,]

intron = as.integer(grepl("-",df$HGVSc))
df=df[intron!=1,]
#remove stop lost and stop gain
stop=as.integer(grepl('Ter',df$HGVSp))

df=df[stop!=1,]

intron = as.integer(grepl("\\+",df$HGVSc))
df=df[intron!=1,]

ns = as.integer(grepl("%",df$HGVSp))
df=df[ns!=1,]

delin=as.integer(grepl("del",df$HGVSc))
df=df[delin!=1,]

stopcodon=as.integer(grepl("\\*",df$HGVSc))
df=df[stopcodon!=1,]

inv=as.integer(grepl("inv*",df$HGVSc))
df=df[inv!=1,]
#table(df$pathogenic)

write.table(df,file=output_file,row.names=FALSE,sep="\t",quote=FALSE,col.names =TRUE)
