################################################################################
#          prs.R      - hypertrophic cardiomyopathy                            # 
################################################################################

# Last updated March 2022.
# Using R Studio Version 1.4.1103, 2009-2021 RStudio, PBC.
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS

# Purpose of script: This script generates polygenic risk scores (PRS) for 
#                    hypertrophic cardiomyopathy for the cohorts of interest.
#                    The input required for this script are the results from the 
#                    extract.variants.by.coordinate.sh script. 
#                    (2) Additionally, this script annotates rows with groups 
#                    status; control, case: monogenic causative variant ABSENT, 
#                    case: monogenic causative variant PRESENT. 
#
#                    The HCM PGS file has effect alleles as reference alleles.
#                    This script tries to retrieve missing variants from VCF files
#                    and assign them homozygous WT GT and the corresponding risk 
#                    effect weight.

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



write.table(rare.disease.analysis.cohort, 
  file = "rare.disease.analysis.cohort.txt", quote = F,
  col.names = T, row.names = T)

rare.disease.analysis.cohort <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/hypertrophic_cardiomyopathy/rare.disease.analysis.cohort.txt",
                                  header = T, stringsAsFactors = F)

####------------ Create a results dataframe from vcf output ----------------####


# Import results txt files into R studio as dataframe
results.grch37.snps <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/hypertrophic_cardiomyopathy/results.grch37.snps.txt",
                               header = T, stringsAsFactors = F)
results.grch38.snps <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/hypertrophic_cardiomyopathy/results.grch38.snps.txt",
                                  header = T, stringsAsFactors = F)


# Prepare results dataframes for merging (match data types)
results.grch37.snps$CHROM <- as.character(results.grch37.snps$CHROM) 


# bind results for each build. 
results.variants.grch38 <- results.grch38.snps
results.variants.grch37 <- results.grch37.snps

# clean up
# rm(results.grch37.snps)
# rm(results.grch38.snps)

#########################################################################################
#########################################################################################

# Note: Polygenic risk scores generated separately for build 38 and 37 ; 
# specifying a chromosomal location in later code as an identifier. 

## Create a genome build 38 version of the pgs scoring file with effect weights.
temp.pgs.grch38.regions <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/hypertrophic_cardiomyopathy/pgs.grch38.regions.txt",
                                      header = F, stringsAsFactors = F) # Import into R
temp.pgs.grch38.regions <- temp.pgs.grch38.regions[, c("V1", "V2")] # extract relevant columns.
pgs000739.grch37 <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/hypertrophic_cardiomyopathy/PGS000739.txt",
                               header = T, stringsAsFactors = F) # Import into R
temp.pgs000739.grch37 <- pgs000739.grch37[, c("rsID", "effect_allele", "other_allele", "effect_weight")] # extract relevant columns.
pgs000739.grch38 <- bind_cols(temp.pgs.grch38.regions, temp.pgs000739.grch37) # bind columns.
pgs000739.grch38 <- plyr::rename(pgs000739.grch38, c("V1" = "chr_name",    # rename() function of plyr
                                                     "V2" = "chr_position"))

# clean up
rm(temp.pgs.grch38.regions)
rm(temp.pgs000739.grch37)


data <- results.variants.grch37
data[data$GT != "1/1" & data$ALT != pgs000739.grch37$effect_allele, ] # keep rows where the value in column GT is not 1/1 and alt

#########################################################################################
########################### GRCh38 #######################################################

homozygousWT <- c("rs1048302", 
                  "rs7556984",
                  "rs2003585",
                  "rs62177303",
                  "rs13061705",
                  "rs4894803",
                  "rs66761011",
                  "rs10052399",
                  "rs3176326",
                  "rs12212795",
                  "rs9320939",
                  "rs7003871",
                  "rs734638",
                  "rs11196085",
                  "rs72840788",
                  "rs1390519",
                  "rs1480036",
                  "rs7301677",
                  "rs41306688",
                  "rs1814880",
                  "rs8033459",
                  "rs7210446",
                  "rs118060942",
                  "rs4799426",
                  "rs117710064",
                  "rs2832230",
                  "rs2070458");
homozygousREF <- c("T",
                   "G",
                   "T",
                   "C",
                   "C",
                   "A",
                   "A",
                   "T",
                   "G",
                   "G",
                   "G",
                   "C",
                   "C",
                   "T",
                   "G",
                   "A",
                   "T",
                   "C",
                   "A",
                   "T",
                   "C",
                   "G",
                   "C",
                   "A",
                   "C",
                   "G",
                   "A")
temp <- results.variants.grch38;
variants <- unique(temp$ID);
samples <- unique(temp$SAMPLE);
for(i in 1:length(homozygousWT)){
  p <- which(temp$ID == homozygousWT[i]); # p= missing rsids.
  if(!length(p)){ # if p has no length (is empty)
    v <- samples;
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
      temp2$ALT <- temp2$REF[1]; # copy ref into alt position.
    }
    temp2$SAMPLE <- v;
    temp2$GT <- "1/1";
    temp <- rbind(temp, temp2);
  }
}
results.variants.grch38 <- temp;

####----------- Generate polygenic risk score for GRCh38 group -------------####


# Note: Polygenic risk scores generated separately for build 38 and 37 ; 
# specifying a chromosomal location in later code as an identifier. 

## Create a genome build 38 version of the pgs scoring file with effect weights.
temp.pgs.grch38.regions <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/hypertrophic_cardiomyopathy/pgs.grch38.regions.txt",
                                       header = F, stringsAsFactors = F) # Import into R
temp.pgs.grch38.regions <- temp.pgs.grch38.regions[, c("V1", "V2")] # extract relevant columns.
pgs000739.grch37 <- read.table("/home/rtaylor3/re_gecip/shared_allGeCIPs/rtaylor3/hypertrophic_cardiomyopathy/PGS000739.txt",
                                      header = T, stringsAsFactors = F) # Import into R
temp.pgs000739.grch37 <- pgs000739.grch37[, c("rsID", "effect_allele", "other_allele", "effect_weight")] # extract relevant columns.
pgs000739.grch38 <- bind_cols(temp.pgs.grch38.regions, temp.pgs000739.grch37) # bind columns.
pgs000739.grch38 <- plyr::rename(pgs000739.grch38, c("V1" = "chr_name",    # rename() function of plyr
                                 "V2" = "chr_position"))

# clean up
rm(temp.pgs.grch38.regions)
rm(temp.pgs000739.grch37)


## Create a numericGT column. GT is defined as character. Take out / and turn 
# each number into numeric character. Then you can add them together.
results.variants.grch38$numericGT <- apply(results.variants.grch38, 1, function(x){
  sum(as.numeric(strsplit(x[9], split = "/")[[1]]));
})

## Annotate each risk allele called with it's effect weight, new column called 'effect'.
# Multiply effect_weight by numericGT.
# In PGS scoring file, effect allele can be the REF (wild type allele). These results will not appear in vcf results.

results.variants.grch38$effect <- 0;
for(i in 1:nrow(pgs000739.grch38)){
  id <- which(results.variants.grch38$CHROM == pgs000739.grch38$chr_name[i] &
                results.variants.grch38$POS == pgs000739.grch38$chr_position[i])
  if(length(id)){
    results.variants.grch38$effect[id] <- results.variants.grch38$numericGT[id]*pgs000739.grch38$effect_weight[i];
  }
}

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



##################################################################################


homozygousWT <- c("rs1048302", 
                  "rs7556984",
                  "rs2003585",
                  "rs62177303",
                  "rs13061705",
                  "rs4894803",
                  "rs66761011",
                  "rs10052399",
                  "rs3176326",
                  "rs12212795",
                  "rs9320939",
                  "rs7003871",
                  "rs734638",
                  "rs11196085",
                  "rs72840788",
                  "rs1390519",
                  "rs1480036",
                  "rs7301677",
                  "rs41306688",
                  "rs1814880",
                  "rs8033459",
                  "rs7210446",
                  "rs118060942",
                  "rs4799426",
                  "rs117710064",
                  "rs2832230",
                  "rs2070458");
homozygousREF <- c("T",
                   "G",
                   "T",
                   "C",
                   "C",
                   "A",
                   "A",
                   "T",
                   "G",
                   "G",
                   "G",
                   "C",
                   "C",
                   "T",
                   "G",
                   "A",
                   "T",
                   "C",
                   "A",
                   "T",
                   "C",
                   "G",
                   "C",
                   "A",
                   "C",
                   "G",
                   "A"
)
temp <- results.variants.grch37;
variants <- unique(temp$ID);
samples <- unique(temp$SAMPLE);
for(i in 1:length(homozygousWT)){
  p <- which(temp$ID == homozygousWT[i]);
  if(!length(p)){
    v <- samples;
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
      temp2$ALT <- temp2$REF[1]; # copy ref into alt position.
    }
    temp2$SAMPLE <- v;
    temp2$GT <- "1/1";
    temp <- rbind(temp, temp2);
  }
}
results.variants.grch37 <- temp;




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
for(i in 1:nrow(pgs000739.grch37)){
  id <- which(results.variants.grch37$CHROM == pgs000739.grch37$chr_name[i] &
                results.variants.grch37$POS == pgs000739.grch37$chr_position[i])
  if(length(id)){
    results.variants.grch37$effect[id] <- results.variants.grch37$numericGT[id]*pgs000739.grch37$effect_weight[i];
  }
}


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
count(collatedeffect, genome_build)
# genome_build   n
# 1       GRCh37  22
# 2       GRCh38 943



################################################################################
# clean up
rm(collatedeffect.grch37)
rm(collatedeffect.grch38)
