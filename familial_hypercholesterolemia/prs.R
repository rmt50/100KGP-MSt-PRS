################################################################################
#            familial hypercholesterolemia  prs.R                              # 
################################################################################

# Last updated March 2022.
# Using R Studio Version 1.4.1103, 2009-2021 RStudio, PBC.
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS

# Purpose of script: This script generates polygenic risk scores (PRS) for 
#                    disease for the familial hypercholesterolemia cohort.
#                    The input required for this script are the results from the 
#                    extract.variants.by.coordinate.sh script. 
#                    (2) Additionally, this script annotates rows with groups 
#                    status; control, case: monogenic causative variant ABSENT, 
#                    case: monogenic causative variant PRESENT. 

#                    Need a bespoke script to incorporate the haplotypes!

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

####------------ Create a results dataframe from vcf output ----------------####


# Import results txt files into R studio as dataframe
results.grch37.snps <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/familial_hypercholesterolemia/results.grch37.snps.txt",
                               header = T, stringsAsFactors = F)
results.grch38.snps <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/familial_hypercholesterolemia/results.grch38.snps.txt",
                                  header = T, stringsAsFactors = F)


# Prepare results dataframes for merging (match data types)
results.grch37.snps$CHROM <- as.character(results.grch37.snps$CHROM) 

results.variants.grch37 <- results.grch37.snps
results.variants.grch38 <- results.grch38.snps


# clean up
rm(results.grch37.snps)
rm(results.grch38.snps)

####----------- Generate polygenic risk score for GRCh38 group -------------####


# Note: Polygenic risk scores generated separately for build 38 and 37 ; 
# specifying a chromosomal location in later code as an identifier. 

## PGS000814 scoring file downloaded from PGS catlog.
## saved and formatted.
## Import into R
pgsfile <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/familial_hypercholesterolemia/PGS000814.txt",
                                      header = T, stringsAsFactors = F) 

#######################################################################################

# code for assigning homozygous reference allele if variant is not called in vcf. 
homozygousWT <- c("rs2479409", "rs4299376", "rs1564348", "rs1800562", "rs3757354", "rs6511720" );
homozygousREF <- c("G", "G", "T", "G", "C", "G") # effect alleles which are reference alleles. 
temp <- results.variants.grch38;
variants <- unique(temp$ID); # list of unique rsIDs
samples <- unique(temp$SAMPLE); # list of unique sampleids in results file. 
for(i in 1:length(homozygousWT)){ # make p vector for all occurences of rsid of of interest
  p <- which(temp$ID == homozygousWT[i]);
  if(!length(p)){
    v <- samples; # put in v vector all samples that do not have the rsid of interest present. 
  }else{
    v <- temp$SAMPLE[p];
    v <- samples[which(!(samples %in% v))];
  }
  if(length(v)){
    if(length(v) == length(samples)){
      temp2 <- temp[rep(1,length(v)),];
      temp2$ALT <- homozygousREF[i];
    }else{
      temp2 <- temp[rep(p[1],length(v)),];
      temp2$ALT <- temp2$REF[1];
    }
    temp2$SAMPLE <- v;
    temp2$GT <- "1/1";
    temp <- rbind(temp, temp2);
  }
}
results.variants.grch38 <- temp;

#####################################################################################
## Commands for applying effect weights for haplotypes.
## Create a dataframe aposamples with assigns haplotype to each sample id. 
temp <- data.frame("p1" = c("T", "T", "T", "T"),
                   "p2" = c("T", "T", "T", "C"),
                   "p3" = c("T", "T", "C", "C"),
                   "p4" = c("T", "C", "T", "C"),
                   "p5" = c("T", "C", "C", "C"),
                   "p6" = c("C", "C", "C", "C"));
apogenotypes <- data.frame("apogenotype" = c("e2;e2", "e2;e3", "e2;e4", "e3;e3",
                                             "e3;e4", "e4;e4"),
                           "string" = apply(temp, 2, function(x)paste(x,collapse  = "")));

# add effect weights to haplotypes. 
apogenotypes$effect_weights <- c(-0.9, -0.4, 0.2, 0, 0.1, 0.2)

apogenotypes$string <- paste(apogenotypes[,2:7])

# Create a dataframe called temp which has all rs429358 and rs7412 calls.
temp <- results.variants.grch38; 
temp <- temp[which(temp$ID == "rs429358" | temp$ID == "rs7412"),];
samples <- unique(temp$SAMPLE)
aposamples <- data.frame("sample" = samples, "genotype" = "")
for(i in 1:length(samples)){
  temp2 <- temp[which(temp$SAMPLE == samples[i]),];
  rs429358 <- c();
  rs7412 <- c();
  for(j in 1:nrow(temp2)){
    if(temp2$ID[j] == "rs429358"){
      if(temp2$GT[j] == "1/1"){
        rs429358 <- c(temp2$ALT[j], temp2$ALT[j]);
      }else if(temp2$GT[j] == "0/1"){
        rs429358 <- c(temp2$REF[j], temp2$ALT[j]);
      }
    }else if(temp2$ID[j] == "rs7412"){
      if(temp2$GT[j] == "1/1"){
        print("This shouldnt happen...check your input data");
        rs7412 <- c(temp2$ALT[j], temp2$ALT[j]);
      }else if(temp2$GT[j] == "0/1"){
        rs7412 <- c(temp2$REF[j], temp2$ALT[j]);
      }
    }
  }
  if(!length(rs7412)){
    rs7412 <- c("C", "C");
  }
  if(!length(rs429358)){
    rs429358 <- c("T", "T");
  }
  if(length(rs7412 > 0) & length(rs429358) > 0){
    aposamples$genotype[i] <- apogenotypes$apogenotype[which(apogenotypes$string == 
                                                               paste(c(rs429358[1], rs7412[2],rs429358[2], rs7412[1]),
                                                                     collapse = ""))];
  }else{
    aposamples$genotype[i] <- NA;
  }
}

## add effect weight to haplotypes. 


####################################################################################

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
for(i in 1:nrow(pgsfile)){
  x <- which(results.variants.grch38$ID == pgsfile$rsID[i] &
               results.variants.grch38$REF == pgsfile$other_allele[i] &
               results.variants.grch38$ALT == pgsfile$effect_allele[i])# value x reflects variables that should match.
  if(length(x)){
    results.variants.grch38$effect[x] <- results.variants.grch38$numericGT[x]*pgsfile$effect_weight[i];
  }else{
    x <- which(results.variants.grch38$ID == pgsfile$rsID[i] &
                 results.variants.grch38$REF != pgsfile$other_allele[i] &
                 results.variants.grch38$ALT == pgsfile$effect_allele[i])
    if(length(x)){
      unmatchedref <- c(unmatchedref, x);
    }
    x <- which(results.variants.grch38$ID == pgsfile$rsID[i] &
                   results.variants.grch38$REF == pgsfile$other_allele[i] &
                   results.variants.grch38$ALT != pgsfile$effect_allele[i])
    if(length(x)){
      unmatchedalt <- c(unmatchedalt, x);
    }
  }
}

## USER CHECK: 
print("Cases where reference allele does not match other allele:");
print(results.variants.grch38[unmatchedref,]);
print("Cases where alternate allele does not match effect allele:");
print(results.variants.grch38[unmatchedalt,]);

count(rare.disease.analysis.cohort, "assembly") 
# assembly freq
# 1   GRCh37    1
# 2   GRCh38  787
#####################################################################################
####################################################################################

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

#################################################################################

# Add on effect weight from haplotypes. 
# Add effect weight column to aposamples and fill.
for(i in 1:nrow(aposamples)){
  aposamples$effect_weight[i] <- apogenotypes$effect_weights[
    which(apogenotypes$apogenotype == aposamples$genotype[i])][1];
}
collatedeffect.grch38_2 <- collatedeffect.grch38

# Add effect weight to summed effect weight in collatedeffect table. 
################ how



####----------- Generate polygenic risk score for GRCh37 group -------------####


# Note: Polygenic risk scores should be generated separately for build 38 and 37,
# as it's be specifying a chromosomal location in later code as an identifier. 

## Create a numericGT column. GT is defined as character. Take out / and turn 
# each number into numeric character. Then you can add them together.
results.variants.grch37$numericGT <- apply(results.variants.grch37, 1, function(x){
  sum(as.numeric(strsplit(x[9], split = "/")[[1]]));
})


## Annotate each risk allele called with it's effect weight, new column called 'effect'.
# Multiply effect_weight by numericGT.
# In PGS scoring file, effect allele can be the REF (wild type allele). These results will not appear in vcf results.
# Hence, created 'unmatchedref_37' and 'unmatchedalt_37' to check for if this happens. 
results.variants.grch37$effect <- 0;
unmatchedref_37 <- c(); # code to allow manual inspection of results where PGS scoring file other allele is not REF allele in participant data. 
unmatchedalt_37 <- c(); # code to allow manual inspection of results where PGS scoring file effect allele is not ALT allele in participant data. 
for(i in 1:nrow(pgsfile)){
  y <- which(results.variants.grch37$POS == pgsfile$chr_position[i] &
               results.variants.grch37$CHROM == pgsfile$chr_name[i] &
               results.variants.grch37$REF == pgsfile$other_allele[i] &
               results.variants.grch37$ALT == pgsfile$effect_allele[i])# value y reflects variables that should match.
  if(length(y)){
    results.variants.grch37$effect[y] <- results.variants.grch37$numericGT[y]*pgsfile$effect_weight[i];
  }else{
    y <- which(results.variants.grch37$POS == pgsfile$chr_position[i] &
                 results.variants.grch37$CHROM == pgsfile$chr_name[i] &
                 results.variants.grch37$REF != pgsfile$other_allele[i] &
                 results.variants.grch37$ALT == pgsfile$effect_allele[i])
    if(length(y)){
      unmatchedref_37 <- c(unmatchedref_37, y);
    }
    y <- which(results.variants.grch37$POS == pgsfile$chr_position[i] &
                 results.variants.grch37$CHROM == pgsfile$chr_name[i] &
                 results.variants.grch37$REF == pgsfile$other_allele[i] &
                 results.variants.grch37$ALT != pgsfile$effect_allele[i])
    if(length(y)){
      unmatchedalt_37 <- c(unmatchedalt_37, y);
    }
  }
}
## USER CHECK: 
print("Cases where reference allele does not match other allele:");
print(results.variants.grch37[unmatchedref_37,]);
print("Cases where alternate allele does not match effect allele:");
print(results.variants.grch37[unmatchedalt_37,]);

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
# assembly freq
# 1   GRCh37    1
# 2   GRCh38  787



################################################################################
# clean up
rm(collatedeffect.grch37)
rm(collatedeffect.grch38)
