################################################################################
#                      extract.variants.by.coordinate.sh                       #               
################################################################################

# Last updated March 2022. 

# Purpose of script: This script subsets VCFs by coordinate and prints out 
#                     tab-delimted files of sample ID, variant information and 
#                     quality.  
#                     !! Script is executed using the login nodes as a portal to
#                     the HPC via a submission as a job: 
#                     <bsub < extract.variants.by.coordinate.sh> !!
#                     Execute in working directory where script is saved.
#                     Login nodes not designed to run large jobs; HPC is used.
#                     Script is subdivided by genome builds in VCF file list and 
#                     chromosome nomenculature.
#                     Script is subdivied with code altered for extraction of 
#                     snps verus indels. 
#                     Tip: To restate job# type <bjob>. To see job progress type 
#                     <bjobs -l [insert job#]> .


## Package information: 
#bash --version echo $
# GNU bash, version 4.4.20(1)-release (x86_64-pc-linux-gnu)
# Copyright (C) 2016 Free Software Foundation, Inc.
# License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>


####------------------------extract pgs alleles  --------------------------####


# STEPS:
# 1: tabix to subset each VCF by the coordinates given. tabix indexes a tab 
#    de-limited file (regions file) i.e. sorts positional data (puts it in 
#    number order). bcftools requires indexed files.
# 2: bcftools norm: to split multi-allelic variants across multiple lines. 
#    Tidies indels. Input for bcftools is the call file (the bcf file).
# 2:'bcftools norm -m -any': -m means split multiallelic sites into biallelic 
#    records, and SNPs and indels should be merged into a single record 
#    (specified by 'any').
# 3: bcftools view: -f PASS to filter for PASS variants; include only sites 
#    which have no filters set.
# 3: For SNPs: filter flag filters for variants with a coverage >10 and genotype 
#    quality >15. 
# 3: For indels: no filter flag as it impedes indel calling.
# 4: bcftools query -f: to print out the sample ID, and variant and quality 
#    information as a text file. 
# 4: bcftolls query: Extracts fields from VCF or BCF files and outputs them in 
#    user-defined format.-f : print. /: use tab seperators. 

#!/bin/bash
#BSUB -q inter
#BSUB -P re_gecip_machine_learning
#BSUB -o /re_gecip/shared_allGeCIPs/rtaylor3/familial_breast_and_or_ovarian_cancer/extract.variants.by.coordinate.%J.out
#BSUB -e /re_gecip/shared_allGeCIPs/rtaylor3/familial_breast_and_or_ovarian_cancer/extract.variants.by.coordinate.%J.err
#BSUB -J extract.variants.by.coordinate.sh
#BSUB -R "rusage[mem=10000] span[hosts=1]"
#BSUB -M 2000
#BSUB -n 2
#BSUB -cwd /re_gecip/shared_allGeCIPs/rtaylor3/familial_breast_and_or_ovarian_cancer


module load bio/BCFtools/1.10.2-GCC-8.3.0


####----------code of build GRCh37 vcf paths and regions file  -------------####


# To extract SNPs:
while read -r vcf; do
tabix -h $vcf -R pgs.grch37.regions.txt | \
bcftools norm -m -any | \
bcftools view -f PASS -i 'MIN(FMT/DP)>10 & MIN(FMT/GQ)>15' | \
bcftools query -f '[%SAMPLE]\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT]\t[%GQ]\t[%DP]\n' >> results.grch37.snps.txt;
done < vcf.file.paths.grch37.txt

sed -i '1s/^/SAMPLE\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tGQ\tDP\n/' results.grch37.snps.txt

# To extract indels only (bcftools view -v indels): 
while read -r vcf; do
tabix -h $vcf -R pgs.grch37.regions.txt | \
bcftools norm -m -any | \
bcftools view -f PASS | \
bcftools view -v indels | \
bcftools query -f '[%SAMPLE]\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT]\t[%GQ]\t[%DP]\n' >> results.grch37.indels.txt;
done < vcf.file.paths.grch37.txt

sed -i '1s/^/SAMPLE\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tGQ\tDP\n/' results.grch37.indels.txt


####---------- code of build GRCh38 vcf paths and regions file -------------####


# To extract SNPs:
while read -r vcf; do
tabix -h $vcf -R pgs.grch38.regions.txt | \
bcftools norm -m -any | \
bcftools view -f PASS -i 'MIN(FMT/DP)>10 & MIN(FMT/GQ)>15' | \
bcftools query -f '[%SAMPLE]\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT]\t[%GQ]\t[%DP]\n' >> results.grch38.snps.txt;
done < vcf.file.paths.grch38.txt

sed -i '1s/^/SAMPLE\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tGQ\tDP\n/' results.grch38.snps.txt

# To extract indels only (bcftools view -v indels): 
while read -r vcf; do
tabix -h $vcf -R pgs.grch38.regions.txt | \
bcftools norm -m -any | \
bcftools view -f PASS | \
bcftools view -v indels | \
bcftools query -f '[%SAMPLE]\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT]\t[%GQ]\t[%DP]\n' >> results.grch38.indels.txt;
done < vcf.file.paths.grch38.txt

sed -i '1s/^/SAMPLE\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tGQ\tDP\n/' results.grch38.indels.txt


################################################################################
