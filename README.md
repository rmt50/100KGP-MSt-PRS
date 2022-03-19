# 100KGP-MSt-PRS
Code for 100KGP MSt project: Assessing the use of polygenic risk scores in the Genomic Medicine Service.

The scripts and files submitted are for use in my MSt Genomic Medicine research project; made available to Genomic Medicine Course Team for assessment purposes. 

(1) Repository contains the folder 'familial_breast_cancer'. This contains scripts which form a pipeline for assessing the polygenic risk for familial breast cancer of certain cohorts within the 100KGP. The pipeline is performed within the Genomics England Research Environment. 

Pipeline step (associated script):
- Create cohort (cohorts.R)
- Create and format files in preparation for extracting risk alleles (vcf.paths.regions.R)
- Extract risk alleles from cohort of interest VCF files (extract.variants.by.coordinates.sh)
- Generate polygenic risk scores for cohort participants (prs.R)
- Format files, add predicted variable information and categorise cohort by groups (predicted.variable.groups.R)
- Perform statistical analysis (statistics.R)

Files supplied in repository (purpose):
- PGS000004.txt (input for script)
- pgs.grch37.regions.txt (input for script)
- pgs.grch38.regions.txt (input for script)
- extract.variants.by.coordinate.851010.out (output of script for examiner interest)
- extract.variants.by.coordinate.851010.err (output of script for examiner interest)


(2) Repository contains the folder 'GSA'. This contains scripts which form a pipeline for assessing variants which are on the GSA. 

Pipeline step (associated script):
- Find risk loci present on GSA, generate files in format suitable for input into pipeline (1) above (gsa.R).

Files supplied in repository (purpose):
- PGS000004.txt (input for script)


(3) Repository contains the folders 'Hypertrophic_Cardiomyopathy' and 'familial_hypercholesterolemia'. 
These contain prs.R scripts that are specific to how to calculate PRS for those diseases.
Special consideration/additional steps needed to be included due to reference alleles being effect alleles and different variants being assigned different weight depending on haplotype. Scripts are not able to run, without the additional files and cohort. Scripts are for examiner to see extra consideration coding steps needed and that this must be considered when calculating PRS.


END
