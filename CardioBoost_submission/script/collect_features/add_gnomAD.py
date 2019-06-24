"""
Script for adding afZscore given a table file 
""" 
import argparse
from collections import defaultdict
import gzip
import pysam
import sys

NEEDED_AF_FIELDS = [ 'gnomAD_AF']
NEEDED_AF_FIELDS_SET = set(NEEDED_AF_FIELDS)
AF_EMPTY_COLUMN_VALUES = ['']*len(NEEDED_AF_FIELDS_SET)
p = argparse.ArgumentParser()
p.add_argument("-i", "--input_table", help="input.tsv", required=True)
p.add_argument("-p", "--af_file", help="icc_missense_gnomAD_AF.txt.gz", required=True)
args = p.parse_args()

counts = defaultdict(int)

def get_af_column_values(af_f, chrom, pos, ref, alt):
    """Retrieves the af value corresponding to the given chrom, pos, ref, alt

    Args:
      af_f: A pysam.TabixFile object corresponding to the gnomAD AF file
      chrom: chromosome (eg. '1')
      pos: the minrepped variant position
      ref: the minrepped ref allele
      ref: the minrepped alt allele

    Return:
      A list of chrom, pos, ref, alt
    """

    if chrom == 'MT':
        return af_EMPTY_COLUMN_VALUES

    counts['total_variants'] += 1

    position_found = False
    af_alt_alleles = []
    for af_f_row in af_f.fetch(chrom, pos-1, pos):
        af_row_fields = af_f_row.split('\t')
        if str(pos) !=  af_row_fields[1]:
            continue
        position_found = True
        af_ref_allele = af_row_fields[2]
        af_alt_allele = af_row_fields[3]
        #print af_row_fields;
        
        if ref == af_ref_allele and alt == af_alt_allele:
            counts['input_variants_with_matching_position_and_matching_allele'] += 1
            break
        af_alt_alleles.append(af_alt_allele)
    else:
        if not position_found:
            counts['input_variants_with_no_matching_position_in_af_file'] += 1
        else:     
#            if ref != af_ref_allele and alt != af_alt_allele:
#                counts['input_variants_with_mismatching_ref_and_alt_allele_in_affile'] += 1   
            if ref != af_ref_allele:
                counts['input_variants_with_mismatching_ref_allele_in_af_file'] += 1
            elif alt != af_alt_allele:
                counts['input_variants_with_mismatching_alt_allele_in_af_file'] += 1
            else:
                counts['input_variants_with_unknown_mismatch'] += 1
            #sys.stderr.write("WARNING: af variant (%s:%s %s>%s) mismatches the input allele (%s:%s %s>%s)\n" % (af_row_fields[0],af_row_fields[1],af_row_fields[3],af_row_fields[4], chrom, pos, ref, alt))
        return AF_EMPTY_COLUMN_VALUES
    
    af_column_values=[af_row_fields[8]];
       
    return af_column_values


af_f = pysam.TabixFile(args.af_file)
input_f = gzip.open(args.input_table) if args.input_table.endswith('.gz') else open(args.input_table)
input_header = next(input_f).rstrip('\n').split('\t')
input_with_af_header = input_header + NEEDED_AF_FIELDS
print("\t".join(input_with_af_header))
for i, input_row in enumerate(input_f):
    input_fields = input_row.rstrip('\n').split('\t')
    input_dict = dict(zip(input_header, input_fields))
    chrom = input_dict['CHROM']
    pos = int(input_dict['POS'])
    ref = input_dict['REF']
    alt = input_dict['ALT']
    af_column_values = get_af_column_values(af_f, chrom, pos, ref, alt)
    
    print("\t".join(input_fields + af_column_values))

for k, v in counts.items():
    sys.stderr.write("%30s: %s\n" % (k, v))

