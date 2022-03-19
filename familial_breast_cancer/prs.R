################################################################################
#                                    prs.R                                     # 
################################################################################

# Last updated March 2022.
# Using R Studio Version 1.4.1103, 2009-2021 RStudio, PBC.
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS

# Purpose of script: This script generates polygenic risk scores (PRS) for 
#                    familial breast cancer for the cohorts of interest.
#                    The input required for this script are the results from the 
#                    extract.variants.by.coordinate.sh script. 
#                    (2) Additionally, this script annotates rows with groups 
#                    status; control, case: monogenic causative variant ABSENT, 
#                    case: monogenic causative variant PRESENT. 

print(sessionInfo(), local = FALSE) # state packages and versions.


####------------------------- Load Packages --------------------------------####


library(tidyverse) # Contains readr, stringr, dplyr packages
library(plyr) # use for re-name function. 

# # Version of packages when last updated/run script.
# Rlabkey_2.8.2  
# httr_1.4.2 
# jsonlite_1.7.3 
# tidyverse_1.3.1
# plyr_1.8.6 
###############################################################################
write.table(rare.disease.analysis.cohort, 
            file = "rare.disease.analysis.cohort.txt", quote = F,
            col.names = T, row.names = T)

rare.disease.analysis.cohort <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/familial_breast_and_or_ovarian_cancer/rare.disease.analysis.cohort.txt",
                                           header = T, stringsAsFactors = F)



####------------ Create a results dataframe from vcf output ----------------####


# Import results txt files into R studio as dataframe
results.grch37.snps <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/familial_breast_and_or_ovarian_cancer/results.grch37.snps.txt",
                               header = T, stringsAsFactors = F)
results.grch37.indels <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/familial_breast_and_or_ovarian_cancer/results.grch37.indels.txt",
                                 header = T, stringsAsFactors = F)
results.grch38.snps <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/familial_breast_and_or_ovarian_cancer/results.grch38.snps.txt",
                                  header = T, stringsAsFactors = F)
results.grch38.indels <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/familial_breast_and_or_ovarian_cancer/results.grch38.indels.txt",
                                    header = T, stringsAsFactors = F)

# Prepare results dataframes for merging (match data types)
results.grch37.snps$DP <- as.character(results.grch37.snps$DP) 
results.grch38.snps$DP <- as.character(results.grch38.snps$DP)
results.grch37.snps$CHROM <- as.character(results.grch37.snps$CHROM) 
results.grch37.indels$CHROM <- as.character(results.grch37.indels$CHROM)

# bind results for each build. 
results.variants.grch37 <- bind_rows(results.grch37.snps, results.grch37.indels)
results.variants.grch38 <- bind_rows(results.grch38.snps, results.grch38.indels)


# clean up
rm(results.grch37.snps)
rm(results.grch38.snps)
rm(results.grch37.indels)
rm(results.grch38.indels)

####----------- Generate polygenic risk score for GRCh38 group -------------####


# Note: Polygenic risk scores generated separately for build 38 and 37 ; 
# specifying a chromosomal location in later code as an identifier. 

## Create a genome build 38 version of the pgs scoring file with effect weights.
temp.pgs.grch38.regions <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/familial_breast_and_or_ovarian_cancer/pgs.grch38.regions.txt",
                                       header = F, stringsAsFactors = F) # Import into R
temp.pgs.grch38.regions <- temp.pgs.grch38.regions[, c("V1", "V2")] # extract relevant columns.
temp.pgs000004.grch37 <- pgs000004.grch37[, c("effect_allele", "other_allele", "effect_weight")] # extract relevant columns.
pgs000004.grch38 <- bind_cols(temp.pgs.grch38.regions, temp.pgs000004.grch37) # bind columns.
pgs000004.grch38 <- plyr::rename(pgs000004.grch38, c("V1" = "chr_name",    # rename() function of plyr
                                 "V2" = "chr_position"))

# clean up
rm(temp.pgs.grch38.regions)
rm(temp.pgs000004.grch37)


## Create a numericGT column. GT is defined as character. Take out / and turn 
# each number into numeric character. Then you can add them together.
results.variants.grch38$numericGT <- apply(results.variants.grch38, 1, function(x){
  sum(as.numeric(strsplit(x[9], split = "/")[[1]]));
})

## Annotate each risk allele called with it's effect weight, new column called 'effect'.
# Multiply effect_weight by numericGT.
# In PGS scoring file, effect allele can be the REF (wild type allele). These results will not appear in vcf results.
# Hence, created 'unmatchedref' and 'unmatchedalt' to check for if this happens. 
results.variants.grch38$effect <- 0;
unmatchedref <- c(); # code to allow manual inspection of results where PGS scoring file other allele is not REF allele in participant data. 
unmatchedalt <- c(); # code to allow manual inspection of results where PGS scoring file effect allele is not ALT allele in participant data.
unmatchedaltANDref <- c(); # code to allow manual inspection of results where PGS scoring file effect allele is not ALT + other is not REF allele in participant data. 
for(i in 1:nrow(pgs000004.grch38)){
  x <- which(results.variants.grch38$POS == pgs000004.grch38$chr_position[i] &
                 results.variants.grch38$CHROM == pgs000004.grch38$chr_name[i] &
               results.variants.grch38$REF == pgs000004.grch38$other_allele[i] &
               results.variants.grch38$ALT == pgs000004.grch38$effect_allele[i]) # makes sure effect applied to effect alelle, not a multialleic base. 
  if(length(x)){
    results.variants.grch38$effect[x] <- results.variants.grch38$numericGT[x]*pgs000004.grch38$effect_weight[i];
  }else{
    x <- which(results.variants.grch38$POS == pgs000004.grch38$chr_position[i] &
                 results.variants.grch38$CHROM == pgs000004.grch38$chr_name[i] &
                 results.variants.grch38$REF != pgs000004.grch38$other_allele[i] &
                 results.variants.grch38$ALT == pgs000004.grch38$effect_allele[i])
    if(length(x)){
      unmatchedref <- c(unmatchedref, x);
    }
    x <- which(results.variants.grch38$POS == pgs000004.grch38$chr_position[i] &
                   results.variants.grch38$CHROM == pgs000004.grch38$chr_name[i] &
                   results.variants.grch38$REF == pgs000004.grch38$other_allele[i] &
                   results.variants.grch38$ALT != pgs000004.grch38$effect_allele[i])
    if(length(x)){
      unmatchedalt <- c(unmatchedalt, x);
    }
    x <- which(results.variants.grch38$POS == pgs000004.grch38$chr_position[i] &
                 results.variants.grch38$CHROM == pgs000004.grch38$chr_name[i] &
                 results.variants.grch38$REF != pgs000004.grch38$other_allele[i] &
                 results.variants.grch38$ALT != pgs000004.grch38$effect_allele[i])
    if(length(x)){
      unmatchedaltANDref <- c(unmatchedaltANDref, x);
    }
  }
}

options(max.print=1000000)
## USER CHECK: 
print("Cases where reference allele does not match other allele:");
unmatchedref.df <- print(results.variants.grch38[unmatchedref,])
print("Cases where alternate allele does not match effect allele:");
unmatchedalt.df <- print(results.variants.grch38[unmatchedalt,]);
print("Cases where alternate/reference allele does not match effect/other allele:");
unmatchedaltANDref.df <- print(results.variants.grch38[unmatchedaltANDref,]);

#############################################################################################
## Tidy: Delete all rows that have 0 in effect column. 
count(results.variants.grch38[which(results.variants.grch38$effect==0), ])
results.variants.grch38 <-subset(results.variants.grch38, effect!="0")
count(results.variants.grch38[which(results.variants.grch38$effect==0), ])
#######################################################################################################################################
## No need to account for homozygous WT as all effect alleles are alternate alelles for PGS0000004.
########################################################################################################


count(rare.disease.analysis.cohort, "assembly") # = 759
# assembly freq
# 1   GRCh37   50
# 2   GRCh38  846
count(rare.disease.analysis.cohort, participant_id, which(rare.disease.analysis.cohort$assembly == "GRCh38")) # = 759
# ## Count how many unique ppids in 38 results.
count(results.variants.grch38, duplicated(results.variants.grch38$SAMPLE)) # = 758
# ################################################################################################################## why is 1 sample dropped???????




################################################################################################################################################
## Create new dataframe named collatedeffect.grch38, with column in it named 
# sample_id' and 'PRS', and 'participant_id'.
collatedeffect.grch38 <- data.frame(
  "sample_id" = unique(results.variants.grch38$SAMPLE), 
  "PRS" = 0, 
  "participant_id" = "");

## Fill in PRS column, which the sum of all the effect scores for each variant, 
# for each participant (named as their LP# from SAMPLE column). 
for(i in 1:nrow(collatedeffect.grch38)){
  collatedeffect.grch38$PRS[i] <- sum(results.variants.grch38$effect[which(
    results.variants.grch38$SAMPLE == collatedeffect.grch38$sample_id[i])]);
  ## Fill in 'participant_id' column in collatedeffect.grch38 dataframe, using input form rare.disease.analysis.cohort, using 'sample_id' (LP#) as the common matching variable as it is in both dataframes.
  collatedeffect.grch38$participant_id[i] <- rare.disease.analysis.cohort$participant_id[which(
    rare.disease.analysis.cohort$plate_key == collatedeffect.grch38$sample_id[i])][1];
}

collatedeffect.grch38$participant_id[duplicated(collatedeffect.grch38$participant_id)]
collatedeffect.grch38$sample_id[duplicated(collatedeffect.grch38$sample_id)]


####----------- Generate polygenic risk score for GRCh37 group -------------####


# Note: Polygenic risk scores should be generated separately for build 38 and 37,
# as it's be specifying a chromosomal location in later code as an identifier. 

## Create a numericGT column. GT is defined as character. Take out / and turn 
# each number into numeric character. Then you can add them together.
results.variants.grch37$numericGT <- apply(results.variants.grch37, 1, function(x){
  sum(as.numeric(strsplit(x[9], split = "/")[[1]]));
})

results.variants.grch37$effect <- 0;
unmatchedref_37 <- c(); # code to allow manual inspection of results where PGS scoring file other allele is not REF allele in participant data. 
unmatchedalt_37 <- c(); # code to allow manual inspection of results where PGS scoring file effect allele is not ALT allele in participant data.
unmatchedaltANDref_37 <- c(); # code to allow manual inspection of results where PGS scoring file effect allele is not ALT + other is not REF allele in participant data. 
for(i in 1:nrow(pgs000004.grch37)){
  y <- which(results.variants.grch37$POS == pgs000004.grch37$chr_position[i] &
               results.variants.grch37$CHROM == pgs000004.grch37$chr_name[i] &
               results.variants.grch37$REF == pgs000004.grch37$other_allele[i] &
               results.variants.grch37$ALT == pgs000004.grch37$effect_allele[i]) # makes sure effect applied to effect alelle, not a multialleic base. 
  if(length(y )){
    results.variants.grch37$effect[y ] <- results.variants.grch37$numericGT[y ]*pgs000004.grch37$effect_weight[i];
  }else{
    y  <- which(results.variants.grch37$POS == pgs000004.grch37$chr_position[i] &
                  results.variants.grch37$CHROM == pgs000004.grch37$chr_name[i] &
                  results.variants.grch37$REF != pgs000004.grch37$other_allele[i] &
                  results.variants.grch37$ALT == pgs000004.grch37$effect_allele[i])
    if(length(y )){
      unmatchedref_37 <- c(unmatchedref_37, y);
    }
    y <- which(results.variants.grch37$POS == pgs000004.grch37$chr_position[i] &
                 results.variants.grch37$CHROM == pgs000004.grch37$chr_name[i] &
                 results.variants.grch37$REF == pgs000004.grch37$other_allele[i] &
                 results.variants.grch37$ALT != pgs000004.grch37$effect_allele[i])
    if(length(y)){
      unmatchedalt_37 <- c(unmatchedalt_37, y);
    }
    y <- which(results.variants.grch37$POS == pgs000004.grch37$chr_position[i] &
                 results.variants.grch37$CHROM == pgs000004.grch37$chr_name[i] &
                 results.variants.grch37$REF != pgs000004.grch37$other_allele[i] &
                 results.variants.grch37$ALT != pgs000004.grch37$effect_allele[i])
    if(length(y)){
      unmatchedaltANDref_37 <- c(unmatchedaltANDref_37, y);
    }
  }
}
## USER CHECK: 
print("Cases where reference allele does not match other allele:");
unmatchedref_37.df <- print(results.variants.grch37[unmatchedref_37,]);
print("Cases where alternate allele does not match effect allele:");
unmatchedalt_37.df <- print(results.variants.grch37[unmatchedalt_37,]);
print("Cases where alternate/reference allele does not match effect/other allele:");
unmatchedaltANDref_37 <- print(results.variants.grch37[unmatchedaltANDref_37,]);

######################################################################################################
# For indels beware of equivalent deletions not getting an effect score assigned when they should. 

###################################################################################################
# Times when effect = 0, is because the position is multialleic, so variant called but reference/alt does not match other/effect. 
# It's good these are assigned 0. 

#############################################################################################
## Tidy: Delete all rows that have 0 in effect column. 
count(results.variants.grch37[which(results.variants.grch37$effect==0), ])
results.variants.grch37 <-subset(results.variants.grch37, effect!="0")
count(results.variants.grch37[which(results.variants.grch37$effect==0), ])
#######################################################################################################################################
## No need to account for homozygous WT as all effect alleles are alternate alelles for PGS0000004.
#######################################################################################################################################


## Create new dataframe named collatedeffect.grch37, with column in it named 
# 'sample_id' and 'PRS', and 'participant_id'.
collatedeffect.grch37 <- data.frame("sample_id" = unique(results.variants.grch37$SAMPLE),
                                    "PRS" = 0,
                                    "participant_id" = "");
## Fill in PRS column, which the sum of all the effect scores for each variant, 
# for each participant (named as their LP# from SAMPLE column). 
for(i in 1:nrow(collatedeffect.grch37)){
  collatedeffect.grch37$PRS[i] <- sum(results.variants.grch37$effect[which(
    results.variants.grch37$SAMPLE == collatedeffect.grch37$sample_id[i])]);
  ## Fill in 'participant_id' column in collatedeffect.grch37 dataframe, 
  # using input form rare.disease.analysis.cohort, using 'sample_id' (LP#) as the common matching variable as it is in both dataframes.
  collatedeffect.grch37$participant_id[i] <- rare.disease.analysis.cohort$participant_id[which(
    rare.disease.analysis.cohort$plate_key == collatedeffect.grch37$sample_id[i])][1];
}




# clean up
rm(pgs000004.grch37)
rm(pgs000004.grch38)
rm(results.variants.grch37)
rm(results.variants.grch38)


####------------------ Merge collatedeffect dataframes ---------------------####

# There may be some participants that have same result in genome build 37 and genome build 38. 
# For these participants, select results from genome build 38 as it is the most 
# updated version of genome; i.e. more instances/chances there can be alignment 
# at PGS allele specific loci.

# collatedeffect.grch38. 758 obs. 3 var.
# collatedeffect.grch37. 33 obs. 3 var. 


collatedeffect <- collatedeffect.grch38;
collatedeffect$genome_build <- "GRCh38";
collatedeffect.grch37$genome_build <- "GRCh37";
collatedeffect <- rbind(collatedeffect,
                        collatedeffect.grch37[which(!(
                          collatedeffect.grch37$sample_id %in% collatedeffect$sample_id)),])

# collatedeffect. 791 obs. 4 var. 
#Number of entries in each status assignment.
count(collatedeffect, "genome_build")
# genome_build freq
# 1       GRCh37   50
# 2       GRCh38  843



################################################################################
# clean up
rm(collatedeffect.grch37)
rm(collatedeffect.grch38)

