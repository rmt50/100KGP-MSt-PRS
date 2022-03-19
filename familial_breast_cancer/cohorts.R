###############################################################################
#                 participant selection     (cohorts.R)                       #
###############################################################################

# Last updated March 2022.
# R version 4.0.3 (2020-10-10)
# R Studio Version 1.4.1103, 2009-2021 RStudio, PBC.

# Purpose of script: This script creates 3 sub-cohorts from the 100KGP data:
#                    cases affected and tested for breast carcinoma panel where 
#                    100KGP results (1) explained phenotype, (2) did not explain 
#                    phenotype, and (3) control non-proband group not affected by 
#                    breast carcinoma. 

## Information of user:
print(sessionInfo(), local = FALSE) # state packages and versions.
####---------------------------  Load Packages --------------------------####

library(jsonlite) # The Rlabkey package also depends on the "jsonlite" package.
library(httr) # The Rlabkey package also depends on the "httr" package.
library(Rlabkey) # LabKey Remote API for R package, which can be obtained via 
# CRAN using the package name "Rlabkey". 
# See https://www.labkey.org/Documentation/21.7/wiki-page.view?name=rAPI for more information.
library(tidyverse) # Contains readr, dplyr, stringr used later.
library(dplyr)
# Version of packages when last updated/run script.
# Rlabkey_2.8.2
# jsonlite_1.7.3
# httr_1.4.2
# tidyverse_1.3.1

####-------------------------------- CASES --------------------------------####

# START with 75,526 Rare Disease participants. 
#Filter for participant ids that are:
# Proband
# Female
# Multiancestry
# Rare Disease programme.
# Not withdrawn
# Has HPO HP:0003002
# Does not have HP:0030075;HP:0030076
# Has Familial breast cancer panel applied. 
# Not null value for age_of_onset for familial breast and or ovarian cancer.
# Particpant ids not under same family-id. 
# Subdivide into causative variant + not. 


participant <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="participant", 
  viewName="", 
  colSelect="participant_id,programme,participant_ethnic_category,participant_phenotypic_sex,programme_consent_status,rare_diseases_family_id,participant_type,duplicated_participant_id", 
  # Various filters to subset rows
  colFilter=makeFilter(
    c("participant_type", "EQUAL", "Proband"),
    c("participant_phenotypic_sex", "EQUAL", "Female"),
    c("programme", "EQUAL", "Rare Diseases"),
    c("participant_id", "NOT_IN", "[insert withdrawn ppids here")), 
  # withdrawn.
  containerFilter=NULL, 
  colNameOpt="rname"
)
# 10327 obs. of 8 var.

count(participant, duplicated(participant$participant_id))
# 10,327 participants.


phenotype <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="rare_diseases_participant_phenotype", 
  viewName="", 
  colSelect="participant_id,hpo_term,hpo_id,laterality,hpo_present", 
  colFilter=makeFilter(
    c("hpo_present", "EQUAL", "Yes"),
    c("hpo_id", "EQUAL", "HP:0003002")), 
  containerFilter=NULL, 
  colNameOpt="rname"
)
#1094 obs. 5 var. 
count(phenotype, duplicated(phenotype$participant_id))
# 1,094 participants.


## Create participants list to filter against.
# participants with HPO terms for in situ cancers.
unwanted_hpos <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="rare_diseases_participant_phenotype", 
  viewName="", 
  colSelect="participant_id,hpo_term,hpo_id,hpo_present", 
  colFilter=makeFilter(
    c("hpo_present", "EQUAL", "Yes"),
    c("hpo_id", "IN", "HP:0030075;HP:0030076")), # 
  containerFilter=NULL, 
  colNameOpt="rname"
)
# 41 obs. 4 var.
count(unwanted_hpos, duplicated(unwanted_hpos$participant_id))
#39 participants

# remove insitu breast cancers.
phenotype_subset <- phenotype[-which(
  phenotype$participant_id %in% unwanted_hpos$participant_id), ]
#1068 obs. 5 var.
count(phenotype_subset, duplicated(phenotype_subset$participant_id))
# 1,068 participants (26 removed) with HPO term of interest. 
# Doesn't need to be 39. Not all pps with HP:0030075;HP:0030076 will also have HP:0003002.

# filter for meeting HPO are participant demographic filters applied ealier. 
participant_filteredforHPO <- merge(
  participant, phenotype_subset, by = "participant_id")
# 516 obs. 12 var. 
count(participant_filteredforHPO, duplicated(participant_filteredforHPO$participant_id))
# 516 participants. 
head(participant_filteredforHPO)

# select rows from panels applied on labkey that had familial breast cancer panel.
panels_applied <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="panels_applied", 
  viewName="", 
  colSelect="participant_id,panel_name,panel_version", 
  colFilter=makeFilter(c("panel_name", "EQUAL", "Familial breast cancer")), 
  containerFilter=NULL, 
  colNameOpt="rname"
)
# 1321 obs. of 3 variables. 
count(panels_applied, duplicated(panels_applied$participant_id))
# 1304 participants. 
panels_applied <- panels_applied %>%  distinct()
# 1304. obs. 3 var. 

# Filter for pps affected with BC who had BC panel applied.
# Important step to check as HPO term presence does not mean that it was main 
# indication/reason for recruitment. 
participant_filteredforPanel <- merge(
  participant_filteredforHPO, panels_applied, by = "participant_id")
#491 obs. 14 variables.
count(participant_filteredforPanel, duplicated(
  participant_filteredforPanel$participant_id))
# 491 participants. 

## Sucessfully filtered for:
# Proband
# Female
# Multiancestry
# Rare Disease programme.
# Not withdrawn
# Has HPO HP:0003002
# Does not have HP:0030075;HP:0030076
# Has Familial breast cancer panel applied. 



## Add disease serverity variable. 
disease <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="rare_diseases_participant_disease", 
  viewName="", 
  colSelect="participant_id,age_of_onset,normalised_specific_disease", 
  colSort="-age_of_onset,-normalised_age_of_onset", 
  colFilter=makeFilter(
    c("age_of_onset", "NOT_IN", "-999;0;999"),
    c("normalised_specific_disease", "EQUAL", "Familial breast and or ovarian cancer")), 
  containerFilter=NULL, 
  colNameOpt="rname"
)
# 808 obs. 3 var. 

## Add 'age_of_onset' column from disease onto the end of the 
# participant_filteredforPanel dataframe, using 'participant_id' as the 
# matching common column variable.
for(i in 1:nrow(participant_filteredforPanel)){
  participant_filteredforPanel$age_of_onset[i] <- disease$age_of_onset[
    which(disease$participant_id == participant_filteredforPanel$participant_id[i])][1];
}

for(i in 1:nrow(participant_filteredforPanel)){
  participant_filteredforPanel$normalised_specific_disease[i] <- disease$normalised_specific_disease[
    which(disease$participant_id == participant_filteredforPanel$participant_id[i])][1];
}
# 491 obs. 16 var. 

count(participant_filteredforPanel, age_of_onset, sort = T)
# 181 N/As. Need to be removed. 

participant_filteredforPanel_agefilter <- drop_na(participant_filteredforPanel, age_of_onset)
# 310 obs. 16 variables.
count(participant_filteredforPanel_agefilter, duplicated(
  participant_filteredforPanel_agefilter$participant_id))
# 310 participants.


## check for relatedness
count(participant_filteredforPanel_agefilter, duplicated(
  participant_filteredforPanel_agefilter$rare_diseases_family_id))
# no duplicates of rare_diseases_family_id. 


## check for duplicate samples (duplicated_participant_id column)
count(participant_filteredforPanel_agefilter, duplicated(
  participant_filteredforPanel_agefilter$duplicated_participant_id))
participant_filteredforPanel_agefilter$duplicated_participant_id[
  duplicated(participant_filteredforPanel_agefilter$duplicated_participant_id)]
# none. 


# ## Categorise cases based on whether a monogenic causative variant has been identified. 
# Select rows, with specified filters, from the rare_diseases_interpreted Labkey table. 
# Data for all interpreted rare disease participants
# filters define cases that have a monogenic causative variant identified. 
# interpretations that has been fully evaluated by a GMC and will appear in the gmc_exit_questionnaire.
# entriy from the gmc_exit_questionnaire =  case_solved_family
## From the guide "Researchers are encouraged to work on the subset of samples 
# that have already passed our internal QC checks; these can be found below for 
# rare disease and cancer genomes, respectively.For Rare Disease genomes, it 
# should be noted that all tiered genomes have passed through Genomics England 
#in-house QCs. 

interpreted <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="rare_disease_interpreted", 
  viewName="", 
  colSelect="participant_id,plate_key,delivery_id,assembly,alignment_file_path,platypus_vcf_path,case_solved_family", 
  colFilter=makeFilter(c("case_solved_family", "IN", "yes;no")), 
  containerFilter=NULL, 
  colNameOpt="rname"
)
# 63411 obs. 7 var. 


# merge tables, filter for pps with vcf that passed internal qc .
cases <- merge(
  participant_filteredforPanel_agefilter, interpreted, by = "participant_id")
# 292 obs. 22 var. 
count(cases, duplicated(
  cases$participant_id))
#292 participants
interpreted$participant_id[duplicated(interpreted$participant_id)]

# make a unique ID column; will use to filter later.
cases$ID <- paste(
  cases$participant_id, 
  cases$plate_key, 
  cases$delivery_id,
  cases$assembly,
  sep = "_");

# Make new column status. 
# Assign case's a causative/no causative monogenic variant status. 
# Assign values to a new column based on values of another column where a condition is satisfied.
cases$status <- "case_monogenic_variant_ABSENT" # default placeholder entry.
cases$status[which((
  cases$case_solved_family == "yes")
  )] <- "case_monogenic_variant_PRESENT"
# 292 obs. 24 var. 

count(cases, status) # observe counts of case type. 

cases_cut <- select(cases, c(
  "participant_id", "plate_key", "ID", "assembly", "age_of_onset", "status"))
# 292 obs. 6 var. 


####------------------------------- CONTROLS -------------------------------####


# Filter for:
# be relative, mother. 
# be female.
# multiancesry
# In rare disease cohort. 
# consenting.
# need to create same ID. 
# not have HP:0003002
# not be under Inherited cancer domain. 

participant_2 <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="participant", 
  viewName="", 
  colSelect="participant_id,programme,participant_ethnic_category,participant_phenotypic_sex,programme_consent_status,rare_diseases_family_id,participant_type,duplicated_participant_id", 
  # Various filters to subset rows
  colFilter=makeFilter(
    c("participant_type", "EQUAL", "Relative"),
    c("participant_phenotypic_sex", "EQUAL", "Female"),
    c("programme", "EQUAL", "Rare Diseases"),
    c("participant_id", "NOT_IN", "[insert withdrawn ppids here]"), # withdrawn.
    c("biological_relationship_to_proband", "EQUAL", "Mother")), 
  containerFilter=NULL, 
  colNameOpt="rname"
)
# 11750 obs. 8 var. 
count(participant_2, duplicated(participant_2$participant_id))
# 11750 participants. 

# create list to filter against.
phenotype_2 <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="rare_diseases_participant_phenotype", 
  viewName="", 
  colSelect="participant_id,hpo_term,hpo_id,hpo_present", 
  colFilter=makeFilter(
    c("hpo_present", "EQUAL", "Yes"),
    c("hpo_id", "EQUAL", "HP:0003002")), 
  containerFilter=NULL, 
  colNameOpt="rname"
)
# 1139 obs. 4 var.
count(phenotype_2, duplicated(phenotype_2$participant_id))
# 1139 participants.

participant_2_filtered <- participant_2[-which(
  participant_2$participant_id %in% phenotype_2$participant_id), ]
# 11673 obs. 8 var
count(participant_2_filtered, duplicated(participant_2_filtered$participant_id))
# 11673 participants. 

## Sucessfully filtered for:
# Relative, mother
# Female
# Rare Disease programme.
# Not withdrawn
# Does not have HPO HP:0003002

## Filter for those without normalised_specific_disease_proband as 
# familial breast and or ovarian cancer. 
# extract entries for columns we'll make ID from. 
interpreted_2 <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v14_2022-01-27", 
  schemaName="lists", 
  queryName="rare_disease_interpreted", 
  viewName="", 
  colSelect="participant_id,plate_key,delivery_id,assembly,normalised_specific_disease_proband,", 
  colFilter=makeFilter(c("normalised_specific_disease_proband", "NOT_EQUAL_OR_MISSING", "Familial breast and or ovarian cancer")), 
  containerFilter=NULL, 
  colNameOpt="rname"
)
# 71324 obs. 5 var. 
count(interpreted_2, duplicated(interpreted_2$participant_id))
# 71315 participants.

# merge
controls <- merge(
  participant_2_filtered, interpreted_2, by = "participant_id")
# 11600 obs. 12 var

count(controls, duplicated(controls$participant_id))
# 11597 participants ids.

controls$participant_id[duplicated(controls$participant_id)]
# participants ids that are dupliacated have entry in build 37+38, #
# or due to different family id. 
controls_unique <- controls %>%  distinct()
controls_unique$participant_id[duplicated(controls_unique$participant_id)]

## check for relatedness
count(controls_unique, duplicated(
  controls_unique$rare_diseases_family_id))
# 2 replicated.
controls_unique$rare_diseases_family_id[duplicated(
  controls_unique$rare_diseases_family_id)]
# inspect. OK, same participant id (same as dups for build in above line)

## check for duplicate samples (duplicated_participant_id column)
count(controls_unique, duplicated(
  controls_unique$duplicated_participant_id))
# 1 entriy that aren't NA. 
controls_unique$duplicated_participant_id[!duplicated(
  controls_unique$duplicated_participant_id)]
# OK, is not in particiany id column. 

# make a unique ID column; will use as a filter later.
controls_unique$ID <- paste(
  controls_unique$participant_id, 
  controls_unique$plate_key, 
  controls_unique$delivery_id,
  controls_unique$assembly,
  sep = "_");

# Make new column status. 
controls_unique$status <- "control" # default placeholder entry.
# 11599 obs. 14 var. 
count(controls_unique, duplicated(controls_unique$participant_id))
# 11597 participants

controls_cut <- select(controls_unique, c(
  "participant_id", "plate_key", "ID", "assembly", "status"))
# 11599 obs. 5 var. 

# Slice dataframe to user defined number of rows to cohort size to suit analysis.
controls500 <- slice(controls_cut, 7001:7500)  # select 500 entries.
count(controls500, duplicated(controls500$participant_id))
# 500 participants (no dup ids.)

####------------------------------- COMBINE-------------------------------####

## Merge
colnames(controls500)
count(controls500, duplicated(controls500$participant_id))
# 500 obs. 5 variables. 
colnames(cases_cut)
count(cases_cut, duplicated(cases_cut$participant_id))
# 292 obs. 6 variables.

rare.disease.analysis.cohort<- bind_rows(cases_cut, controls500, id = NULL)
count(rare.disease.analysis.cohort.2, duplicated(rare.disease.analysis.cohort$participant_id))
# 792 obs. 6 variables.
# 792 participants


# Information for user.
#Number of entries in each status assignment.
count(rare.disease.analysis.cohort, "status")
# status freq
# 1  case_monogenic_variant_ABSENT  379
# 2 case_monogenic_variant_PRESENT   17
# 3                        control  500




count(rare.disease.analysis.cohort, duplicated(rare.disease.analysis.cohort$participant_id))

################################################################################

# clean up
rm(cases, cases_cut, controls, controls_cut, controls_unique, controls500, 
   disease, interpreted, interpreted_2, panels_applied, participant, 
   participant_2, participant_2_filtered, participant_filteredforHPO,
   participant_filteredforPanel, participant_filteredforPanel_agefilter, 
   phenotype, phenotype_2, phenotype_subset, unwanted_hpos)



