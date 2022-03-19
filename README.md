# 100KGP-MSt-PRS
Code for 100KGP MSt project: Assessing the use of polygenic risk scores in the Genomic Medicine Service.

The scripts and files submitted are for use in my MSt Genomic Medicine research project.
Repository is set to private, available to Genomic Medicine Course Team for assessment purposes only. 
The scripts form a pipeline for assessing the polygenic risk of certain diseases for cohorts of interest within the 100KGP.
The pipeline is performed within the Genomics England Research Environment. 

Pipeline step (associated script):
- Create cohort (cohorts.R)
- Create and format files in preparation for extracting risk alleles (vcf.paths.regions.R)
- Extract risk alleles from cohort of interest VCF files (extract.variants.by.coordinates.sh)
- Generate polygenic risk scores for cohort participants (prs.R)
- Format files, add predicted variable information and categorise cohort by groups (predicted.variable.groups.R)
- Perform statistical analysis (statistics.R)

Files supplied in repository:
- PGS000004.txt
- pgs.grch37.regions.txt
- pgs.grch38.regions.txt
