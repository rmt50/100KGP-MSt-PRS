################################################################################
#                             gsa.R                                            #
################################################################################

# Last updated March 2022.

# Purpose of this script: This script finds which 313 BC polygenic risk associated loci
#                         are on the GSAv3 array.
#                         It outputs a tab de-limited files of common loci in build 
#                         37 and 38, some with effect weight/alt/ref info.
#                         Output files made are called:
#                         sharedVariants_arrayAndSequencing_grch38.txt
#                         and sharedVariants_arrayAndSequencing_grch37.txt
#                         and pgs.regions.grch37.txt and pgs.regions.grch38.txt


#### --------------------------- Load packages --------------------------- ####


library(tidyverse)


#### -------------------------- Find common loci -------------------------- ####


## USER ACTION REQUIRED: Download 'Infinium Global Screening Array v3.0 Manifest File (CSV Format – GRCh37)'
# from https://support.illumina.com/downloads/infinium-global-screening-array-v3-0-product-files.html. 
## USER ACTION REQUIRED: Delete header (first 7 lines) and save as GSA-24v3-0_A1_noHeading.txt.

## Read GSA-24v3-0_A1_noHeading.txt into R.
gsafile <- read.csv("GSA-24v3-0_A1_noHeading.txt",header = T, stringsAsFactors = F)
nrow(gsafile) # Check number of rows
head(gsafile) # Visualise file. 
# Columns of note in gdsfile:
# Chr: Chromosome containing the SNP.
# MapInfo: Chromosomal coordinates of the SNP.

## Read PGS catalog scoring file into R. 
pgsfile <- read.table("PGS000004.txt", header = T, stringsAsFactors = T);
head(pgsfile)

## Makes a unique ID by combining Chromosome Chromosomal coordinates information
# for GSA file and the PGS file.
gsafile$ID <- paste(gsafile$Chr, gsafile$MapInfo, sep = "_");
pgsfile$ID <- paste(pgsfile$chr_name, pgsfile$chr_position, sep = "_");


## Cross-compare the unique IDs and saves the shared variants.
# Note: both files are in the genome build GRCh37. 
commonvars <- pgsfile[which(pgsfile$ID %in% gsafile$ID),];
write.table(commonvars, "sharedVariants_arrayAndSequencing.grch37.txt", col.names = T,
            row.names = F, quote = F)

## Make a regions file in build 37
write.table(commonvars[,c("chr_name","chr_position", "chr_position")], "pgs.regions.grch37.txt",
            col.names = F, row.names = F, quote = F) 

## Make a regions file in build 37 compatible with liftover website (add chr prefix)
temp <- commonvars; 
temp$chr <- paste("chr", temp$chr_name, sep = "");
write.table(temp[,c(ncol(temp),2,2)], "sharedVariants_arrayAndSequencing_trimmed.txt",
            col.names = F, row.names = F, quote = F) # create input for LiftOver.

## USER ACTION REQUIRED: Use white listed page https://genome.ucsc.edu/cgi-bin/hgLiftOver
# Convert sharedVariants_arrayAndSequencing_trimmed.txt from build 37 to 38. 
# save output as pgs.regions.grch38.txt 

# read into R the LiftOver output.
temp <- read.table("pgs.regions.grch38.txt", header = F, stringsAsFactors = F); 
commonvars38 <- commonvars;
commonvars$chr_position <- temp$V2;
write.table(commonvars, "sharedVariants_arrayAndSequencing_grch38.txt",
            col.names = T, row.names = F, quote = F)


###### END
