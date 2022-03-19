###########################################################################
#                       predicted.variable.groups.R                       # 
###########################################################################

# Last updated March 2022.
# Using R Studio Version 1.4.1103, 2009-2021 RStudio, PBC.

# Purpose of script: (1) This script adds values for the predicted variable age
#                        onset, to the dataframe.
#                    (2) Additionally, this script annotates rows with groups 
#                        status; control, case: monogenic causative variant ABSENT
#                       , case: monogenic causative variant PRESENT. 
#                    (3) The input required for this script is collatedeffect df.

####---------------------------------- Load packages --------------------------------####

library(jsonlite)
library(httr)
library(Rlabkey) 
library(tidyverse)


####-------------- Add disease severity information (age of onset)----------####

## Add 'age_of_onset' column from rare.disease.analysis.cohort onto the end 
# of the collatedeffect dataframe, using 'participant_id' as the matching common 
# column variable.
for(i in 1:nrow(collatedeffect)){
  collatedeffect$age_of_onset[i] <- rare.disease.analysis.cohort$age_of_onset[
    which(rare.disease.analysis.cohort$participant_id == collatedeffect$participant_id[i])][1];
}

count(collatedeffect, "age_of_onset")


####--------------------Add group status information ----------------------####

colnames(collatedeffect)
colnames(rare.disease.analysis.cohort)

## Annotates rows with group status; control, case: monogenic causative variant 
# ABSENT, case: monogenic causative variant PRESENT.
# These annotations were originally assigned in the rare.disease.analysis.cohort 
# dataframe, so annotation can be carried over. 
# Sample id / platekey used as the matching common column variable.
for(i in 1:nrow(collatedeffect)){
  collatedeffect$status[i] <- rare.disease.analysis.cohort$status[which(rare.disease.analysis.cohort$plate_key == collatedeffect$sample_id[i])[1]];
}

count(collatedeffect, "status")
# status freq
# status freq
# 1  case_monogenic_variant_ABSENT  376
# 2 case_monogenic_variant_PRESENT   17
# 3                        control  500


count(collatedeffect, "genome_build")
# genome_build freq
# 1       GRCh37   50
# 2       GRCh38  843

###################################################################################################

