import gzip
from ngs_parsers import VCF
from collections import OrderedDict

def setup_likelihood_file(fout, vcf):

    l = vcf.empty_vcf_line.keys()[9:]
    header = [item for sublist in zip(l,l,l) for item in sublist]
    header = ['I','SID'] + header

    fout = open(fout, 'w')
    fout.write('\t'.join(header) + "\n")
    fout.close()

def float_2_decimal(value):
    return '{:f}'.format(value)

phred_2_pvalue = lambda x: 10**((-1*x)/10)

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

fin_path = "/Users/testudines/DATA/CAP_MAR.uni_geno.snps.indels.recallibrated.vcf.gz"
fin = gzip.open(fin_path,'rb')

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


# GATK: SNP Qual is is -10 * log(1-p) that a REF/ALT polymorphism exists at this site given sequencing data.
# DBM: SNP quality is used to compute the probability that the SNP is a false positive SNP, by formula 10^(-quality / 30).

fin.close()

vcf = VCF.VCF(input=fin_path)

likelihoods_fout = 'test.liks'

setup_likelihood_file(likelihoods_fout, vcf)

likelihoods_fout = open(likelihoods_fout,'a')
snps_fout = open('test.snps','w')

for c, vcf_line_dict in enumerate(vcf.vcf_file_iterator()):
    #print vcf_line_dict
    chrm =  vcf_line_dict['CHROM']
    chrm_id = chrm_2_id_dict[chrm]
    # SNP file:
    snp_file_line = ('rs{}'.format(c), chrm_id, vcf_line_dict['POS'],
                      vcf_line_dict['QUAL'], vcf_line_dict['REF'],
                      vcf_line_dict['ALT'])

    snps_fout.write("\t".join(map(str, snp_file_line)) + "\n")


    triplets = convert_PLs_to_likelihoods(vcf_line_dict)
    line = [item for sublist in triplets for item in sublist] # flatten triplets
    line = map(str, line)
    line = ['NA','rs{}'.format(c)] + line
    likelihoods_fout.write('\t'.join(line) + '\n')

snps_fout.close()
likelihoods_fout.close()