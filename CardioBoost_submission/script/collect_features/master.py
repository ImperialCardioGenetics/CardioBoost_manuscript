import os
import sys

#path of vep command https://www.ensembl.org/info/docs/tools/vep/index.html
VEPPATH='path/to/vep'

#R script to extract the variant effect annotation on canonical transcripts
EXTRACT_VEP_CANON = "./extract_vep_canon.R"

DISEASESET = ['arm','cm']
ARM_GENE_INFO = "../../data/arrhythmia/arm_gene.txt"
CM_GENE_INFO = "../../data/cardiomyopathy/cm_gene.txt"

#use annovar annotation script http://annovar.openbioinformatics.org/en/latest/user-guide/download/
ANNOVAR = "path/to/annovar/table_annovar.pl"
ANNOVAR_DB = "path/to/annovar/humandb"
CLEAN_ANNOVAR = "./clean_annovar.R"

#http://exac.broadinstitute.org/downloads
ADD_EXAC = "./add_exac_fields.py"
EXAC_FILE = "path/to/ExAC.r0.3.1.sites.vep.normalized.vcf.gz"

#add parazscore https://zenodo.org/record/817898
ADD_PARAZSCORE = "./add_para_zscore.py"
PARAZSCORE_FILE = "path/to/hg19.all.para_zscore.tsv.gz"

#add MPC and misbadness score https://www.biorxiv.org/content/early/2017/06/12/148353
ADD_MPC = "./add_MPC_constraint.py"
MPC_FILE = "path/to/fordist_constraint_official_mpc_values.txt.gz"

MERGE_ANNOVAR_EXAC = "./merge_annovar_exac.py"

ADD_gnomAD = "./add_gnomAD.py"
gnomAD_FILE="path/to/gnomAD_AF_nonpass.txt.gz"


ADD_MSA = "./add_msa.py"

#path of the folder containing the MSA file (both nucleotide and AA)
ARM_GENE_MSA = "../../data/arrhythmia/MSA"
CM_GENE_MSA = "../../data/cardiomyopathy/MSA/"

if __name__ == '__main__':
    disease = sys.argv[1]
    if disease not in DISEASESET:
        sys.exit("The input disease is not in the supported disease set.")
    input_pathogenic_info = sys.argv[2]
    ## The input columns should be in this form: CHROM POS REF ALT pathogenic

    input_file_vcf = sys.argv[3]


    #Option to skip the vep step
    if len(sys.argv) > 4:
        skip_vep = bool(sys.argv[4])
    else:
        skip_vep = False

    if skip_vep == True:
        prefix = input_file_vcf.replace("_vep_canon.vcf","")
    else:
        prefix = input_file_vcf.replace(".vcf","")


    if skip_vep is False:
    #run VEP
        print("Step 1: running VEP")
        vep_command = VEPPATH + " -i " + input_file_vcf + \
                  " --cache --port 3337 --hgvs --symbol --minimal --tab --force_overwrite -o " + prefix +"_vep.tsv"
        print vep_command
        os.system(vep_command)
    #extract the canonical transcripts from vep annotation
        print "Step 2: extracting the canonical transcripts"
        if disease == "arm":
            gene_info = ARM_GENE_INFO
        elif disease == "cm":
            gene_info = CM_GENE_INFO
        extract_command = "Rscript " + EXTRACT_VEP_CANON + " --pathogenic "+ input_pathogenic_info + \
                      " --vep "+ prefix +"_vep.tsv" + " --transcript " + gene_info + \
                      " --out " + prefix + "_vep_canon.txt"
        print extract_command
        os.system(extract_command)

    #Get the annovar annotation
    print "Step 3: running Annovar"
    annovar_command = ANNOVAR + " " + input_file_vcf + " -vcfinput " + ANNOVAR_DB + " -protocol dbnsfp33a,mcap,revel,parazscore -operation f,f,f,f -build hg19 -nastring ."
    print annovar_command
    os.system(annovar_command)
    clean_command = "Rscript " + CLEAN_ANNOVAR + " " + input_file_vcf+".hg19_multianno.txt " + prefix + "_annovar.txt"
    print clean_command
    os.system(clean_command)

    #Join with ExAC
    print "Step 4: getting AF from ExAC"
    exac_command = "python2 -u " + ADD_EXAC + " -i " + prefix + "_vep_canon.txt" + " -e " + EXAC_FILE + " > "+ prefix + "_vep_canon_exac.txt"
    print(exac_command)
    os.system(exac_command)
    os.system("ex -s +'bufdo!v/\S/d' -cxa "+ prefix + "_vep_canon_exac.txt")

    print "Step 5: Merging Annovar and ExAC file"
    merge_command = "python2 " + MERGE_ANNOVAR_EXAC+" "+prefix+"_annovar.txt " + prefix + "_vep_canon_exac.txt " + prefix + "_vep_canon_exac_annovar.txt"
    print(merge_command)
    os.system(merge_command)


    print "Step 6: Adding gnomAD AF"
    gnomAD_command = "python2 " + ADD_gnomAD + " -i " + prefix + "_vep_canon_exac_annovar.txt" + " -p " + gnomAD_FILE + " > " + prefix + "_annotated_1.txt"
    print(gnomAD_command)
    os.system(gnomAD_command)

    print "Step 7: Adding MPC"
    mpc_command = "python2 " + ADD_MPC + " -i " + prefix + "_annotated_1.txt" + " -p " + MPC_FILE + " > " + prefix + "_annotated_2.txt"
    print(mpc_command)
    os.system(mpc_command)


    print "Step 8: Adding features from Orthology MSA"
    if disease == "arm":
        gene_msa = ARM_GENE_MSA
    elif disease == "cm":
        gene_msa = CM_GENE_MSA
    msa_command = "python2 " + ADD_MSA + " -d "+disease + " -i " + prefix + "_annotated_2.txt" + " -m " + gene_msa + " -o " + prefix + "_annotated.txt"
    print(msa_command)
    os.system(msa_command)
