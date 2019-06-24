"""
Script for adding paraZscore given a table file 
""" 
import argparse
from collections import defaultdict
import gzip
import pysam
import sys

NEEDED_PARA_FIELDS = [ 'PARASCOREzSCORE']
NEEDED_PARA_FIELDS_SET = set(NEEDED_PARA_FIELDS)
PARA_EMPTY_COLUMN_VALUES = ['']*len(NEEDED_PARA_FIELDS_SET)
p = argparse.ArgumentParser()
p.add_argument("-i", "--input_table", help="input.tsv", required=True)
p.add_argument("-p", "--para_file", help="hg19.all.para_zscore.tsv.gz", required=True)
args = p.parse_args()

counts = defaultdict(int)

def get_para_column_values(para_f, chrom, pos, ref, alt):
    """Retrieves the paraZscore value corresponding to the given chrom, pos, ref, alt

    Args:
      para_f: A pysam.TabixFile object corresponding to the para score file
      chrom: chromosome (eg. '1')
      pos: the minrepped variant position
      ref: the minrepped ref allele
      ref: the minrepped alt allele

    Return:
      A list of chrom, pos, ref, alt
    """

    if chrom == 'MT':
        return PARA_EMPTY_COLUMN_VALUES

    counts['total_variants'] += 1

    position_found = False
    para_alt_alleles = []
    for para_f_row in para_f.fetch(chrom, pos-1, pos):
        para_row_fields = para_f_row.split('\t')
        if str(pos) !=  para_row_fields[1]:
            continue
        position_found = True
        para_ref_allele = para_row_fields[3]
        para_alt_allele = para_row_fields[4]
        
        #print para_row_fields;
        
        if ref == para_ref_allele and alt == para_alt_allele:
            counts['input_variants_with_matching_position_and_matching_allele'] += 1
            break
        para_alt_alleles.append(para_alt_allele)
    else:
        if not position_found:
            counts['input_variants_with_no_matching_position_in_parafile'] += 1
        else:     
#            if ref != para_ref_allele and alt != para_alt_allele:
#                counts['input_variants_with_mismatching_ref_and_alt_allele_in_parafile'] += 1   
            if ref != para_ref_allele:
                counts['input_variants_with_mismatching_ref_allele_in_parafile'] += 1
            elif alt != para_alt_allele:
                counts['input_variants_with_mismatching_alt_allele_in_parafile'] += 1
            else:
                counts['input_variants_with_unknown_mismatch'] += 1
            #sys.stderr.write("WARNING: Para variant (%s:%s %s>%s) mismatches the input allele (%s:%s %s>%s)\n" % (para_row_fields[0],para_row_fields[1],para_row_fields[3],para_row_fields[4], chrom, pos, ref, alt))
        return PARA_EMPTY_COLUMN_VALUES
    
    para_column_values=[para_row_fields[13]];
       
    return para_column_values


para_f = pysam.TabixFile(args.para_file)
input_f = gzip.open(args.input_table) if args.input_table.endswith('.gz') else open(args.input_table)
input_header = next(input_f).rstrip('\n').split('\t')
input_with_para_header = input_header + NEEDED_PARA_FIELDS
print("\t".join(input_with_para_header))
for i, input_row in enumerate(input_f):
    input_fields = input_row.rstrip('\n').split('\t')
    input_dict = dict(zip(input_header, input_fields))
    chrom = input_dict['CHROM']
    pos = int(input_dict['POS'])
    ref = input_dict['REF']
    alt = input_dict['ALT']
    para_column_values = get_para_column_values(para_f, chrom, pos, ref, alt)
    
    print("\t".join(input_fields + para_column_values))

for k, v in counts.items():
    sys.stderr.write("%30s: %s\n" % (k, v))

