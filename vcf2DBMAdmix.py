#!/usr/bin/env python
# encoding: utf-8


"""
./vcf2DBMAdmix.py \
-i /Users/testudines/Data/CAP_MAR.uni_geno.snps.indels.recallibrated.vcf.gz \
-l CAP_MAR.uni_geno.snps.indels.recallibrated.lik \
-s CAP_MAR.uni_geno.snps.indels.recallibrated.snps \
-r 5


/Users/testudines/Source/dbm/./dbm \
/Users/testudines/Code/vcf_conversion_scripts/CAP_MAR.uni_geno.snps.indels.recallibrated

qsub -V -N dbm -pe omp 2 -l h_rt=4:00:00 -b y \


"""

import gzip
import argparse
from ngs_parsers import VCF
from collections import OrderedDict

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()

    parser.add_argument('-i',
                        '--input',
                        required=True,
                        help='Path to VCF file.')

    parser.add_argument('-l',
                        '--likelihoods-file',
                        help='Path to output file that contains \
                        SNV likelihoods.')

    parser.add_argument('-s',
                        '--SNPs-file',
                        help='Path to output file that contains \
                        SNP positions and ids')

    parser.add_argument('-r',
                        '--region',
                        default=None,
                        type=str,
                        help='chrm:start-stop')

    args = parser.parse_args()

    if args.region is not None:
        if len(args.region.split(":")) == 2:
            chrm = [args.region.split(":")[0]]
            start_stop = [int(item) for item in args.region.split(":")[1].split("-")]
            args.region = chrm + start_stop

    else:
        args.region = [args.region]

    return args

def setup_likelihood_file(fout, vcf):

    l = vcf.empty_vcf_line.keys()[9:]
    header = [item for sublist in zip(l,l,l) for item in sublist]
    header = ['I','SID'] + header

    fout = open(fout, 'w')
    fout.write('\t'.join(header) + "\n")
    fout.close()

def float_2_decimal(value):
    return '{:f}'.format(value)

def phred_2_pvalue(x):
    return 10**((-1*x)/10)

def convert_PLs_to_likelihoods(vcf_line_dict):

    triplets = []

    for i in vcf_line_dict.items()[9:]:

        if i[-1] is None:
            value = (1, 1, 1)  # Not sure if this is the best way to define missing values

        else:
            triplet = i[-1]['PL']
            value = tuple(map(phred_2_pvalue, triplet))

            triplets.append(value)

    return triplets


def convert_chr_names_to_integer_ids(fin):

    chrm_2_id_dict = OrderedDict()
    chrm_2_id_counter = 1

    for c, l in enumerate(fin):
        l = l.strip()

        if l.startswith("##contig"):

            Id = l.split(',')[0].split('=')[-1]

            # Skip unincorporated contigs / scaffolds
            if Id.startswith('A')  or Id.startswith('G'):
                continue

            chrm_2_id_dict[Id] = 'chr{}'.format(chrm_2_id_counter)
            chrm_2_id_counter +=1

        if l.startswith('#CHROM'):
            break

    fin.close()
    return chrm_2_id_dict

if __name__ == '__main__':

    args = get_args()

    vcf = VCF.VCF(input=args.input)
    fin = vcf.__open_vcf__()
    chrm_2_id_dict = convert_chr_names_to_integer_ids(fin)


    setup_likelihood_file(args.likelihoods_file, vcf)
    likelihoods_fout = open(args.likelihoods_file,'a')
    snps_fout = open(args.SNPs_file,'w')

    # Do slicing
    if args.region[0] is None:
        data = vcf.vcf_file_iterator()
    else:
        data = vcf.vcf_slice_iterator(vcf.input, args.region)

    # Iterate over vcf file/slice
    for c, vcf_line_dict in enumerate(data):

        chrm =  vcf_line_dict['CHROM']
        chrm_id = chrm_2_id_dict[chrm]

        # GATK: SNP Qual is -10 * log(1-p) that a REF/ALT polymorphism exists at this site given sequencing data.
        # DBM: SNP quality is used to compute the probability that the SNP is a false positive SNP, by formula 10^(-quality / 30).

        # Skip multi allelic sites
        if len(vcf_line_dict['ALT']) > 1:
            continue

        snp_file_line = ('rs{}'.format(c), chrm_id, vcf_line_dict['POS'],
                          vcf_line_dict['QUAL'] * 3, vcf_line_dict['REF'],
                          vcf_line_dict['ALT'])

        snps_fout.write("\t".join(map(str, snp_file_line)) + "\n")

        triplets = convert_PLs_to_likelihoods(vcf_line_dict)
        line = [item for sublist in triplets for item in sublist] # flatten triplets
        line = map(str, line)
        line = ['NA','rs{}'.format(c)] + line
        likelihoods_fout.write('\t'.join(line) + '\n')

    snps_fout.close()
    likelihoods_fout.close()