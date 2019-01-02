"""
Script for adding mpcZscore given a table file 
""" 
import argparse
from collections import defaultdict
import gzip
import pysam
import sys

NEEDED_mpc_FIELDS = [ 'mis_badness',  'MPC']
NEEDED_mpc_FIELDS_SET = set(NEEDED_mpc_FIELDS)
mpc_EMPTY_COLUMN_VALUES = ['']*len(NEEDED_mpc_FIELDS_SET)
p = argparse.ArgumentParser()
p.add_argument("-i", "--input_table", help="input.tsv", required=True)
p.add_argument("-p", "--mpc_file", help="fordist_constraint_official_mpc_values.txt.gz", required=True)
args = p.parse_args()

counts = defaultdict(int)

def get_mpc_column_values(mpc_f, chrom, pos, ref, alt):
    """Retrieves the mpcZscore value corresponding to the given chrom, pos, ref, alt

    Args:
      mpc_f: A pysam.TabixFile object corresponding to the mpc score file
      chrom: chromosome (eg. '1')
      pos: the minrepped variant position
      ref: the minrepped ref allele
      ref: the minrepped alt allele

    Return:
      A list of chrom, pos, ref, alt
    """

    if chrom == 'MT':
        return mpc_EMPTY_COLUMN_VALUES

    counts['total_variants'] += 1

    position_found = False
    mpc_alt_alleles = []
    for mpc_f_row in mpc_f.fetch(chrom, pos-1, pos):
        mpc_row_fields = mpc_f_row.split('\t')
        if str(pos) !=  mpc_row_fields[1]:
            continue
        position_found = True
        mpc_ref_allele = mpc_row_fields[2]
        mpc_alt_allele = mpc_row_fields[3]
        
        #print mpc_row_fields;
        
        if ref == mpc_ref_allele and alt == mpc_alt_allele:
            counts['input_variants_with_matching_position_and_matching_allele'] += 1
            break
        mpc_alt_alleles.append(mpc_alt_allele)
    else:
        if not position_found:
            counts['input_variants_with_no_matching_position_in_mpcfile'] += 1
        else:     
#            if ref != mpc_ref_allele and alt != mpc_alt_allele:
#                counts['input_variants_with_mismatching_ref_and_alt_allele_in_mpcfile'] += 1   
            if ref != mpc_ref_allele:
                counts['input_variants_with_mismatching_ref_allele_in_MPCfile'] += 1
            elif alt != mpc_alt_allele:
                counts['input_variants_with_mismatching_alt_allele_in_MPCfile'] += 1
            else:
                counts['input_variants_with_unknown_mismatch'] += 1
            #sys.stderr.write("WARNING: mpc variant (%s:%s %s>%s) mismatches the input allele (%s:%s %s>%s)\n" % (mpc_row_fields[0],mpc_row_fields[1],mpc_row_fields[3],mpc_row_fields[4], chrom, pos, ref, alt))
        return mpc_EMPTY_COLUMN_VALUES
    
    mpc_column_values=[mpc_row_fields[16],mpc_row_fields[18]];
       
    return mpc_column_values


mpc_f = pysam.TabixFile(args.mpc_file)
input_f = gzip.open(args.input_table) if args.input_table.endswith('.gz') else open(args.input_table)
input_header = next(input_f).rstrip('\n').split('\t')
input_with_mpc_header = input_header + NEEDED_mpc_FIELDS
print("\t".join(input_with_mpc_header))
for i, input_row in enumerate(input_f):
    input_fields = input_row.rstrip('\n').split('\t')
    input_dict = dict(zip(input_header, input_fields))
    chrom = input_dict['CHROM']
    pos = int(input_dict['POS'])
    ref = input_dict['REF']
    alt = input_dict['ALT']
    mpc_column_values = get_mpc_column_values(mpc_f, chrom, pos, ref, alt)
    
    print("\t".join(input_fields + mpc_column_values))

for k, v in counts.items():
    sys.stderr.write("%30s: %s\n" % (k, v))

