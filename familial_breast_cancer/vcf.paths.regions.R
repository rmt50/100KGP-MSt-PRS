################################################################################
#                         vcf.paths.regions.R                                  #         
################################################################################

# Last updated February 2022.
# Using R Studio Version 1.4.1103, 2009-2021 RStudio, PBC.

# Purpose of script: This script creates the input files required for the 
#                    extract.variants.by.coordinate script.
#                    These files are a list of full-paths to VCF files of 
#                    interest (for GRCh37 and GRCh38) and 
#                    a tab-delimited list of chromosomal coordinates for the 
#                    polygenic risk allele positions (for GRCh37 and GRCh38).
#                    See to the USER ACTION REQUIRED points below prior to 
#                    running script.


####-------------------------- Load Packages -------------------------------####


library(Rlabkey) # LabKey Remote API for R package, which can be obtained via CRAN using the package name "Rlabkey". See https://www.labkey.org/Documentation/21.7/wiki-page.view?name=rAPI for more information.
library(tidyverse) # Contains readr
library(jsonlite) # The Rlabkey package also depends on the "jsonlite" package.
library(httr)


####--------- Create list of full-paths to VCF files for cohort ------------####


# VCF genomic data files are as they have been delivered to GeL by sequencing 
# provider (Illumina). 
# These have all passed an initial QC check based on sequencing quality and coverage. 

# Select rows, with specified filter for vcf file paths, from the 
# genome_file_paths_and_types Labkey table into a data frame called 
# 'temp.genome.file.paths.and.types'.
temp.genome.file.paths.and.types <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="genome_file_paths_and_types", 
  viewName="", 
  colSelect="participant_id,lab_sample_id,platekey,delivery_id,delivery_date,delivery_version,genome_build,type,file_path,filename,file_sub_type,file_type", 
  colFilter=makeFilter(c("file_sub_type", "EQUAL", "Standard VCF")), 
  containerFilter=NULL, 
  colNameOpt="rname"
)
colnames(temp.genome.file.paths.and.types)
# make a unique ID column; will use as a filter later.
temp.genome.file.paths.and.types$ID <- paste(
  temp.genome.file.paths.and.types$participant_id, 
  temp.genome.file.paths.and.types$platekey, 
  temp.genome.file.paths.and.types$delivery_id,
  temp.genome.file.paths.and.types$genome_build,
  sep = "_");

head(temp.genome.file.paths.and.types)


## Add 'file_path' column from disease onto the end of the 
# rare.disease.analaysis.cohort dataframe, using 'ID' as the 
# matching common column variable.
for(i in 1:nrow(rare.disease.analysis.cohort)){
  rare.disease.analysis.cohort$file_path[i] <- temp.genome.file.paths.and.types$file_path[
    which(temp.genome.file.paths.and.types$ID == rare.disease.analysis.cohort$ID[i])][1];
}

# clean up
rm(temp.genome.file.paths.and.types)

# Create vcf pathway list file for each genome build. Requirement for 
# extract_variants script, where input of genome builds in VCF file list cannot 
# be mixed.
# note: Approximately 10% of the genomic data are aligned against the reference 
# genome version GRCh37 and the remaining majority (90%) against version GRCh38.
# Information for user on # of entries in each build in cohort:
count(rare.disease.analysis.cohort, assembly)

write.table(rare.disease.analysis.cohort$file_path[which(
  rare.disease.analysis.cohort$assembly == "GRCh37")], 
  file = "vcf.file.paths.grch37.txt", quote = F,
            col.names = F, row.names = F)
write.table(rare.disease.analysis.cohort$file_path[which(
  rare.disease.analysis.cohort$assembly == "GRCh38")], 
  file = "vcf.file.paths.grch38.txt", quote = F,
            col.names = F, row.names = F)

####--- Create list of chromosomal coordinates for the risk allele loci ----####


# USER ACTION REQUIRED: Download scoring file from PGS catalog 
# webpage: https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000004/ScoringFiles/).
# USER ACTION REQUIRED: Save in working directory as a tab de-limited file (.txt).
# Note: PGS000004.txt has chromosomal positions in GRCh37.

# Import PGS scoring file (PGS000004.txt) into R as dataframe.
pgs000004.grch37 <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/familial_breast_and_or_ovarian_cancer/PGS000004.txt",
                                  header = T, stringsAsFactors = F)

# Create grch37 regions txt file by extracting selected columns in pgs000004.grch37 to create list of regions in compatible input format for extract.variants.by.coordinate.sh
write.table(pgs000004.grch37[,c("chr_name","chr_position", "chr_position")], file="pgs.grch37.regions.txt", quote = F,
          col.names = F, row.names = F)

## Create grch38 regions txt file:
temp <- read.table("home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/familial_breast_and_or_ovarian_cancer/pgs.grch37.regions.txt",
                   header = F, stringsAsFactors = F)
temp$chr <- paste("chr", temp$V1, sep = ""); # liftover required chr prefix to chromosome number. 
write.table(temp[,c(ncol(temp),2,2)], "liftoverinput.txt", quote = F,
          col.names = F, row.names = F, quote = F))

# USER ACTION REQUIRED: import liftoverinput.txt into white listed 
# web page http://www.genome.ucsc.edu/cgi-bin/hgLiftOver and liftover from 
# GRCh37 to GRCh38. 
# USER ACTION REQUIRED: Save converted file into working directory as 
# pgs.grch38.regions.txt.
