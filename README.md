vcf_conversion_scripts
======================

Repo to hold vcf_processing_scripts

### Installation:

Install dependancies:

    brew install samtools
    brew install tabix
    pip install pysam

Clone repo and update submodule:

    $ git clone https://github.com/ngcrawford/vcf_conversion_scripts.git
    $ cd vcf_conversion_scripts/
    $ git submodule update --init

### Instructions for vcf2DBMAdmix.py

Bgzip and index vcf file:

    bgizp my.vcf
    tabix -p vcf my.vcf.gz

Create the `.lik` and `.snps` files:

    cd vcf_conversion_scripts/
    ./vcf2DBMAdmix.py \
    -i path/2/my.vcf.gz \
    -l path/2/my.lik \
    -s path/2/my.snps

Run DBM haplotype caller:

    ./dbm my -like


