"""
Script for generating a new variant table with the ExAC fields below added to each input variant that's in ExAC
"""
import argparse
from collections import defaultdict
import pysam
import sys

NEEDED_EXAC_FIELDS = ['AC_Adj','AN_Adj']

NEEDED_EXAC_FIELDS_SET = set(NEEDED_EXAC_FIELDS)
EXAC_EMPTY_COLUMN_VALUES = ['']*len(NEEDED_EXAC_FIELDS_SET)

p = argparse.ArgumentParser()
p.add_argument("-i", "--variant-table", help="variant .tsv", required=True)
p.add_argument("-e", "--exac-sites-vcf", help="ExAC sites VCF", required=True)
args = p.parse_args()

counts = defaultdict(int)

def get_exac_column_values(exac_f, chrom, pos, ref, alt):
    """Retrieves the ExAC vcf row corresponding to the given chrom, pos, ref, alt, and extracts the column values listed in NEEDED_EXAC_FIELDS

    Args:
      exac_f: A pysam.TabixFile object corresponding to the ExAC vcf
      chrom: chromosome (eg. '1')
      pos: the minrepped variant position
      ref: the minrepped variant ref allele
      alt: the minrepped variant alt allele

    Return:
      A list of chrom, pos, ref, alt
    """

    if chrom == 'MT':
        return EXAC_EMPTY_COLUMN_VALUES

    counts['total_variants'] += 1

    # retrieve ExAC variant - pysam.fetch(..) sometimes returns more than 1 vcf record, so need to filter here
    position_found = False
    exac_alt_alleles = []
    pos_1 = pos-1;
    for exac_vcf_row in exac_f.fetch(str(chrom), pos_1, pos):
        exac_row_fields = exac_vcf_row.split('\t')
        if str(pos) !=  exac_row_fields[1]:
            continue
        position_found = True
        exac_ref_allele = exac_row_fields[3]
        exac_alt_allele = exac_row_fields[4]
        if "," in exac_alt_allele:
            raise Exception("Found multiallelic variant: %s. Expecting an ExAC VCF that has been decomposed / normalized with vt." % "-".join(exac_vcf_row_fields[0:5]))

        if ref == exac_ref_allele and alt == exac_alt_allele:
            counts['variants_with_matching_position_and_matching_allele'] += 1
            break
        exac_alt_alleles.append(exac_alt_allele)
    else:
        if not position_found:
            counts['variants_with_no_matching_position_in_exac'] += 1
        else:
            if len(ref) + len(alt) + len(exac_ref_allele) + len(exac_alt_allele) > 4:
                counts['indel_with_no_matching_allele_in_exac'] += 1
            elif ref != exac_ref_allele and alt != exac_alt_allele:
                counts['snp_with_mismatching_ref_and_alt_allele_in_exac'] += 1
            elif ref != exac_ref_allele:
                counts['snp_with_mismatching_ref_allele_in_exac'] += 1
            elif alt != exac_alt_allele:
                counts['snp_with_mismatching_alt_allele_in_exac'] += 1
            else:
                counts['snp_with_unknown_mismatch'] += 1

            sys.stderr.write("WARNING: ExAC variant %s:%s (http://exac.broadinstitute.org/variant/%s-%s-%s-%s) - ExAC alleles (%s:%s %s>%s) mismatch the variant allele (%s:%s %s>%s)\n" % (chrom, pos, chrom, pos, exac_row_fields[3], exac_row_fields[4], chrom, pos, exac_ref_allele, ",".join(exac_alt_alleles), chrom, pos, ref, alt))

        return EXAC_EMPTY_COLUMN_VALUES

    filter_value = exac_row_fields[6]
    info_fields = [('Filter', filter_value)] + [tuple(kv.split('=')) for kv in exac_row_fields[7].split(';')]
    info_fields = filter(lambda kv: kv[0] in NEEDED_EXAC_FIELDS_SET, info_fields)
    info_fields = dict(info_fields)
    exac_column_values = [info_fields[k] for k in NEEDED_EXAC_FIELDS]

    # check that the variant alt allele matches (one of the) ExAC alt allele(s)
    #if len(alt_alleles) > 1:
    #    # select the AC/AN numbers corresponding to the specific alt allele
    #    alt_allele_index = alt_alleles.index(alt)
    #    exac_column_values = map(lambda x: x.split(",")[alt_allele_index] if "," in x else x, exac_column_values)

    return exac_column_values


exac_f = pysam.TabixFile(args.exac_sites_vcf)
variant_f = open(args.variant_table)
variant_header = next(variant_f).rstrip('\n').split('\t')
variant_with_exac_header = variant_header + NEEDED_EXAC_FIELDS
print("\t".join( variant_with_exac_header))
for i, variant_row in enumerate(variant_f):
    variant_fields = variant_row.rstrip('\n').split('\t')
    variant_dict = dict(zip(variant_header, variant_fields))


    chrom = variant_dict['chrom']

    pos = int(variant_dict['pos'])
    ref = variant_dict['ref']
    alt = variant_dict['alt']

    exac_column_values = get_exac_column_values(exac_f, chrom, pos, ref, alt)

    print("\t".join(variant_fields + exac_column_values))

for k, v in counts.items():
    sys.stderr.write("%30s: %s\n" % (k, v))
