import gzip
from collections import OrderedDict

fin = gzip.open("/Users/testudines/DATA/CAP_MAR.uni_geno.snps.indels.recallibrated.vcf.gz",'rb')

chrm_2_id_dict = OrderedDict()
chrm_2_id_counter = 1

header_line = None

for c, l in enumerate(fin):
    l = l.strip()

    if l.startswith("##contig"):

        Id = l.split(',')[0].split('=')[-1]
        chrm_2_id_dict[Id] = 'chr{}'.format(chrm_2_id_counter)
        chrm_2_id_counter +=1

    if l.startswith('#CHROM'):
        header_line = l.strip().split('\t')
        continue

    if header_line is not None:
        print l
        break
