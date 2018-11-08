import os
import sys

VEPPATH='/Users/xzhang13/Documents/vep/ensembl-vep/vep'
EXTRACT_VEP_CANON = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/Brugada/script/extract_vep_canon.R"

DISEASESET = ['arm','cm']
#BRS_GENE_INFO = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/Brugada/gene_list/brs_icc_mutations_gene.txt"
ARM_GENE_INFO = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/LQTS/gene_list/icc_mutations/lqts_icc_mutations_gene.txt"
CM_GENE_INFO = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/CM/gene_list/cm_icc_mutations_gene.txt"

ANNOVAR = "/Users/xzhang13/Documents/NGStools/annovar/table_annovar.pl"
ANNOVAR_DB = "/Users/xzhang13/Documents/NGStools/annovar/humandb"
CLEAN_ANNOVAR = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/Brugada/script/clean_annovar.R"

ADD_EXAC = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/LQTS/Oxford/add_exac_fields.py"
EXAC_FILE = "/Users/xzhang13/Documents/git/data/exac/ExAC.r0.3.1.sites.vep.normalized.vcf.gz"

ADD_PARAZSCORE = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/LQTS/script/add_para_zscore.py"
PARAZSCORE_FILE = "/Users/xzhang13/Desktop/OneDrive/database/para/hg19.all.para_zscore.tsv.gz"

ADD_MPC = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/LQTS/script/add_MPC_constraint.py"
MPC_FILE = "/Users/xzhang13/Desktop/OneDrive/database/regional_constraint/fordist_constraint_official_mpc_values.txt.gz"

MERGE_ANNOVAR_EXAC = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/Brugada/script/merge_annovar_exac.py"

ADD_gnomAD = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/script/add_feature/add_gnomAD.py"
gnomAD_FILE="/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/Database/gnomAD/icc_missense_gnomAD_AF_nonpass.txt.gz"


#BRS_ADD_MSA = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/script/add_msa/add_new_feature_brs.py"
#ARM_ADD_MSA = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/script/add_msa/add_new_feature_arm.py"
ADD_MSA = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/script/add_msa/add_msa.py"

#BRS_GENE_MSA = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/Brugada/msa/"
ARM_GENE_MSA = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/LQTS/MSA/"
CM_GENE_MSA = "/Users/xzhang13/Desktop/OneDrive/OneDrive\ -\ Imperial\ College\ London/Project1/CM/MSA/"

if __name__ == '__main__':
    disease = sys.argv[1]
    if disease not in DISEASESET:
        sys.exit("The input disease is not in the supported disease set.")
    input_pathogenic_info = sys.argv[2]
    ##TODO have a script to check the input format
    ## CHROM POS REF ALT pathogenic

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
    #os.remove(prefix+"_vep.tsv")



    #Get the annovar annotation
    print "Step 3: running Annovar"
    annovar_command = ANNOVAR + " " + input_file_vcf + " -vcfinput " + ANNOVAR_DB + " -protocol dbnsfp33a,mcap,revel,parazscore -operation f,f,f,f -build hg19 -nastring ."
    print annovar_command
    os.system(annovar_command)
    clean_command = "Rscript " + CLEAN_ANNOVAR + " " + input_file_vcf+".hg19_multianno.txt " + prefix + "_annovar.txt"
    print clean_command
    os.system(clean_command)
    #os.remove(input_file_vcf+".hg19_multianno.txt")

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
    #os.remove(prefix+"_annovar.txt")
    #os.remove(prefix+"_vep_canon_exac.txt")

    print "Step 6: Adding gnomAD AF"
    gnomAD_command = "python2 " + ADD_gnomAD + " -i " + prefix + "_vep_canon_exac_annovar.txt" + " -p " + gnomAD_FILE + " > " + prefix + "_annotated_1.txt"
    print(gnomAD_command)
    os.system(gnomAD_command)

    print "Step 7: Adding MPC"
    #parazscore_command = "python2 " + ADD_PARAZSCORE + " -i " + prefix + "_vep_canon_exac_annovar.txt" + " -p " + PARAZSCORE_FILE + " > " + prefix+"_annotated_1.txt"
    #print parazscore_command
    #os.system(parazscore_command)
    mpc_command = "python2 " + ADD_MPC + " -i " + prefix + "_annotated_1.txt" + " -p " + MPC_FILE + " > " + prefix + "_annotated_2.txt"
    print(mpc_command)
    os.system(mpc_command)
    #os.remove(prefix+"_annotated_1.txt")


    print "Step 8: Adding features from Orthology MSA"
    if disease == "arm":
        gene_msa = ARM_GENE_MSA
    elif disease == "cm":
        gene_msa = CM_GENE_MSA
    msa_command = "python2 " + ADD_MSA + " -d "+disease + " -i " + prefix + "_annotated_2.txt" + " -m " + gene_msa + " -o " + prefix + "_annotated.txt"
    print(msa_command)
    os.system(msa_command)
    #os.remove(prefix + "_annotated_2.txt")

