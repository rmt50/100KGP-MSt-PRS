################################################################################
#                                  statistics.R                                # 
################################################################################

# Last updated March 2022.
# Using R Studio Version 1.4.1103, 2009-2021 RStudio, PBC.

# Purpose of script: This script performs statistical analysis.
#                    The input required for this script is the collatedeffect df.

####------------------------- Load Packages --------------------------------####

library(ggplot2)
library(rstatix)
library(tidyverse) 
library(RColorBrewer)
install.packages("gridExtra")
library("gridExtra")
install.packages("ggResidpanel")
library(ggResidpanel)

write.table(collatedeffect, file = "collatedeffect.txt", quote = F,
            col.names = T, row.names = F)

####----------------------- Summarise / Visualise data ---------------------####
summary(collatedeffect)

# get summary statistics for age of onset and PRS.
summary.by.status <- collatedeffect %>% 
  group_by(status) %>% 
  get_summary_stats(type = "common")
pdf("summary.by.status_gridExtra.pdf") # Export PDF
grid.table(summary.by.status)
dev.off()

## Quality check
# get summary statistics for depth from vcf outputs. 
# using snp specific results tables, not indel as no DP values.
temp.snp.results <- bind_rows(results.grch37.snps, results.grch38.snps)
temp.snp.results$DP <- as.numeric(temp.snp.results$DP) 
summary.depth <- summary(temp.snp.results$DP)
pdf("summary.depth.pdf")      
grid.table(summary.depth)
dev.off()


####----------- Compare age of onset between case status groups ------------####

# get common summary stats for the age_of_onset
collatedeffect[collatedeffect$status != "control",] %>% 
  select(-PRS) %>% 
  group_by(status) %>% 
  get_summary_stats(type = "common")

# visualise in box plot:
collatedeffect[collatedeffect$status != "control",] %>% 
  ggplot(aes(x = status, y = age_of_onset, fill = status)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) + # hiding outliers for confidentiality.
  scale_fill_brewer(palette="Set3") +
  theme_bw()+
  ggtitle("Box plot for status groups of age.of.onset")


#Assumption check: Histogram to assess normality for case_monogenic_variant_ABSENT
collatedeffect[collatedeffect$status == "case_monogenic_variant_ABSENT",] %>% 
  ggplot(aes(x = age_of_onset, colour = status, fill = status)) +
  geom_histogram(bins = 15, colour = "black", show.legend = FALSE) +
  facet_grid(status ~ ., scales = "free") + 
  theme_bw()+
  scale_fill_brewer(palette="Set3") +
  coord_cartesian(ylim = c(7, 70))+ # starting y axis at 7 here shows 5 on graph; for confidentiality. 
  ggtitle("Histograms for status groups of age.of.onset")

#Assumption check: Histogram to assess normality not possible for case_monogenic_variant_PRESENT (can't export as <5 freq for age groups (each yr))
collatedeffect[collatedeffect$status == "case_monogenic_variant_PRESENT",] %>% 
  ggplot(aes(x = age_of_onset, colour = status, fill = status)) +
  geom_histogram(bins = 10, colour = "black", show.legend = FALSE) +
  facet_grid(status ~ ., scales = "free") + 
  theme_bw()+
  scale_fill_brewer(palette="Set3") +
  ggtitle("Histograms for status groups of age.of.onset")   #### CAN'T EXPORT THIS GRAPH

#Assumption check: Construct a Q-Q Plot of the quantiles of the data against the quantiles of a normal distribution:
collatedeffect[collatedeffect$status != "control",] %>% 
  group_by(status) %>%
  ggplot(aes(sample = age_of_onset)) +
  stat_qq() +
  stat_qq_line(colour = "red") +
  facet_wrap(facets = vars(status)
  ) 


# Assumption check: perform a Shapiro-Wilk test on both groups separately.
Shapiro.Wilk_1 <- collatedeffect[collatedeffect$status != "control",] %>% 
  select(-PRS) %>% 
  group_by(status) %>% 
  shapiro_test(age_of_onset)
pdf("Shapiro.Wilk_1_gridExtra.pdf")      
grid.table(Shapiro.Wilk_1)
dev.off()


# Assumption check: Test for equity of variance.
levene_test_1 <- collatedeffect[collatedeffect$status != "control", ] %>% 
  levene_test(age_of_onset ~ status)
pdf("levene_test_1.pdf")      
grid.table(levene_test_1)
dev.off()

###  Bartlett test of homogeneity of variances
temp <- collatedeffect[collatedeffect$status != "control", ]
bartlett.test_1 <- bartlett.test(age_of_onset ~ status,
                                 data = temp)
pdf("bartlett_test_1.pdf")     
grid.table(bartlett_test_1)
dev.off()



# Implement two-sample, two-tailed, t-test:
two_sample.two_tailed.t_test_1 <- collatedeffect[collatedeffect$status != "control",] %>% 
  t_test(age_of_onset ~ status,
         alternative = "two.sided",
         var.equal = TRUE)
pdf("two_sample.two_tailed.t_test_1.pdf")       
grid.table(two_sample.two_tailed.t_test_1)
dev.off()


# If not normally distributed: Implement a two-tailed, Mann-Whitney U test:
wilcox_test_1 <- collatedeffect[collatedeffect$status != "control",] %>% 
  wilcox_test(age_of_onset ~ status,
              alternative = "two.sided")
pdf("wilcox_test_1.pdf")      
grid.table(wilcox_test_1)
dev.off()


#summary data again to get mean values:
collatedeffect[collatedeffect$status != "control",] %>%  
  select(-PRS) %>% 
  group_by(status) %>% 
  get_summary_stats(type = "common")


# ####------------ Compare PRS between status groups----------------####

## ANOVA - multiple samples of continuous data (three or more groups).
# used to find out if the samples came from parent distributions with the same mean. 

# get common summary stats
collatedeffect %>%
  select(-age_of_onset) %>%
  group_by(status) %>%
  get_summary_stats(type = "common")


# visualise in box plot:
collatedeffect %>%
  ggplot(aes(x = status, y = PRS, fill = status)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette="Set3") +
  theme_bw()+
  ggtitle("Box plot for status groups of PRS")
# export of image, adjust size, save

# #Assumption check: histogram of PRS for all status groups.
collatedeffect %>%
ggplot(aes(x = PRS, colour = status, fill = status)) +
geom_histogram(colour = "black", bins = 15, show.legend = FALSE) +
facet_grid(status ~ ., scales = "free") +
theme_bw()+
scale_fill_brewer(palette="Set3") +
ggtitle("Histograms for status groups of PRS")


# Draw overlaying histogram of PRS for all status groups.
ggplot(collatedeffect, aes(x = PRS, colour = status, fill = status)) +
  geom_histogram(colour = "darkgrey", position = "identity", alpha = 0.5, bins = 15) +
  theme_bw()+
  scale_fill_brewer(palette="Set3")+
  ggtitle("Overlay of histograms for status groups of PRS")

#Assumption check: Construct a Q-Q Plot of the quantiles of the data against the quantiles of a normal distribution:
collatedeffect %>%
  group_by(status) %>%
  ggplot(aes(sample = PRS)) +
  stat_qq() +
  stat_qq_line(colour = "red") +
  facet_wrap(facets = vars(status)
  )


#  Assumption check: perform a Shapiro-Wilk test
Shapiro.Wilk_2 <- collatedeffect %>%
  select(-age_of_onset) %>%
  group_by(status) %>%
  shapiro_test(PRS)
pdf("Shapiro.Wilk_2_gridExtra.pdf")
grid.table(Shapiro.Wilk_2)
dev.off()

#Assumption check:create a linear model, extract the residuals and check their normality:
lm_collatedeffect<- lm(PRS ~ status,
                       data = collatedeffect)

# extract the residuals
resid_collatedeffect <- residuals(lm_collatedeffect)

# perform Shapiro-Wilk test on residuals
shapiro_test_3 <- resid_collatedeffect %>% 
  shapiro_test()
pdf("shapiro_test_3.pdf")
grid.table(shapiro_test_3)
dev.off()


###  Bartlett test of homogeneity of variances
bartlett.test_2 <- bartlett.test(PRS ~ status,
              data = collatedeffect)
pdf("bartlett_test_2.pdf")
grid.table(bartlett_test_2)
dev.off()

lm_collatedeffect %>% 
  resid_panel()

# Implement test: Perform an ANOVA test on the data:
anova.lm_collatedeffect <- anova(lm_collatedeffect)
pdf("anova.lm_collatedeffect.pdf")
grid.table(anova.lm_collatedeffect)
dev.off()

## perform Tukey’s range test (or any other post-hoc tests) if the preceding 
# ANOVA test showed that there was a significant difference between the groups 
#Post-hoc testing (Tukey’s rank test)
collatedeffect_tukey <- lm_collatedeffect %>% 
  tukey_hsd()
pdf("collatedeffect_tukey.pdf")
grid.table(collatedeffect_tukey)
dev.off()


#summary data again to get mean values:
collatedeffect %>%
  select(-age_of_onset) %>%
  group_by(status) %>%
  get_summary_stats(type = "common")


###  Kruskal-Wallis is analogous to ANOVA but used when normality not met. 
# Assumption check: Test for equity of variance.
levene_test_2 <- collatedeffect %>%
  levene_test(PRS ~ status)
pdf("levene_test_2.pdf")
grid.table(levene_test_2)
dev.off()

##warning about group coerced to factor. There is no need to worry about this - 
# Levene’s test needs to compare different groups and because aggression is 
# encoded as a numeric value, it converts it to a categorical one before running the test.

# implement Kruskal-Wallis test
collatedeffect.kruskal_test <- collatedeffect %>% 
  kruskal_test(PRS ~ status)
pdf("collatedeffect.kruskal_test.pdf")
grid.table(collatedeffect.kruskal_test)
dev.off()

collatedeffect.kruskal_test_2 <- kruskal.test(PRS ~ status, data = collatedeffect)
pdf("collatedeffect.kruskal_test_2.pdf")
grid.table(collatedeffect.kruskal_test_2)
dev.off()


###
# equivalent of Tukey’s range test for non-normal data is Dunn’s test.
# perform Dunn's test
collatedeffect.dunn_test <- collatedeffect %>% 
  dunn_test(PRS ~ status)
pdf("collatedeffect.dunn_test.pdf")
grid.table(collatedeffect.dunn_test)
dev.off()



###-To test if PRS explains variable penetrance in cases with monogenic var--###


# Make new table for only entries with causative varaint present status. 
causative <- data.frame(filter(collatedeffect, status == "case_monogenic_variant_PRESENT"))

# Summarise and visualise.
# create scatterplot of the data:
causative %>% 
  ggplot(aes(x = PRS, y = age_of_onset)) +
  geom_point()

lm_1 <- lm(age_of_onset ~ PRS,
           data = causative)

# Next, we can create diagnostic plots for the model:
lm_1 %>% 
  resid_panel(plots = c("resid", "qq", "ls", "cookd"),
              smoother = TRUE)
#The diagnostic plots show residuals in four different ways:
# Residual Plot: Residuals vs Fitted. Used to check the linear relationship assumptions. A horizontal line, without distinct patterns is an indication for a linear relationship, what is good.
# Q-Q Plot: Used to examine whether the residuals are normally distributed. It’s good if residuals points follow the straight dashed line.
# Location-Scale Plot: Used to check the homogeneity of variance of the residuals (homoscedasticity). Horizontal line with equally spread points is a good indication of homoscedasticity. 
# Residuals vs Leverage. Used to identify influential cases, that is extreme values that might influence the regression results when included or excluded from the analysis. 
## (1) **Linearity of the data**: The linearity assumption can be checked by inspecting the Residuals vs Fitted plot (1st plot):
plot(lm_1,1) 
## (2) **Homogeneity of variance**: The homogeneity of variance assumption can be checked by inspecting the 'Location-Scale' plot (3rd plot):
plot(lm_1,3) #Closer look at 'Location-Scale' individually.
model2 <- lm(log(age_of_onset) ~ PRS, data = causative) # use a log or square root transformation of the outcome variable (y)
plot(model2, 3)
## (3) **Normality of residuals**: The normality assumption can be checked by inspecting the Q-Q plot (2nd plot).
plot(lm_1,2)
## (4) **Outliers and high levarage points**: Outliers and high leverage points can be identified by inspecting the 'Residuals vs Leverage plot' (secret plot not on 4 plot view, it is the 5th plot)
plot(lm_1,5, id.n=26) 
plot(lm_1,4, id.n=26)
## (5 cont) Also, The Residuals vs Leverage plot can help us to find influential observations if any. On this plot, outlying values are generally located at the upper right corner or at the lower right corner. Those spots are the places where data points can be influential against a regression line.
plot(lm_1,5, id.n=26)

### Implement Test
lm_1


anova_lm1 <- anova(lm_1)
pdf("anova_lm1.pdf")
grid.table(anova_lm1)
dev.off()


#Plotting the regression line
#It can be very helpful to plot the regression line with the original data to see how far the data are from the predicted linear values. We can do this with:
# plot the data
causative %>% 
  ggplot(aes(x = PRS, y = age_of_onset)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("PRS v age_of_onset for cases with a monogenic causative variant present")

causative %>% 
  ggplot(aes(x = PRS, y = age_of_onset)) +
  coord_cartesian(ylim = c(26, 52))+ 
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("PRS v age_of_onset for cases with a monogenic causative variant present") # NO POINTS CAN EXPORT.


###-To test if PRS explains variable penetrance in cases withOUT monogenic var--####################


# Make new table for only entries with causative varaint present status. 
non_causative <- data.frame(filter(collatedeffect, status == "case_monogenic_variant_ABSENT"))

# Summarise and visualise.
# create scatterplot of the data:
non_causative %>% 
  ggplot(aes(x = PRS, y = age_of_onset)) +
  geom_point()

lm_2 <- lm(age_of_onset ~ PRS,
           data = non_causative)

# Next, we can create diagnostic plots for the model:
lm_2 %>% 
  resid_panel(plots = c("resid", "qq", "ls", "cookd"),
              smoother = TRUE)

plot(lm_2,1) 
plot(lm_2,3) 
model3 <- lm(log(age_of_onset) ~ PRS, data = non_causative) 
plot(model3, 3)

plot(lm_2,2)
plot(lm_2,5, id.n=26) 
plot(lm_2,4, id.n=26)
plot(lm_1,5, id.n=26)

### Implement Test
lm_2

anova_lm2 <- anova(lm_2)
pdf("anova_lm2.pdf")
grid.table(anova_lm2)
dev.off()

#Plotting the regression line
non_causative %>% 
  ggplot(aes(x = PRS, y = age_of_onset)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("PRS v age_of_onset for cases with monogenic causative variant absent")

causative %>% 
  ggplot(aes(x = PRS, y = age_of_onset)) +
  coord_cartesian(xlim = c(-1.4, 0.75), ylim = c(17, 73))+ 
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("PRS v age_of_onset for cases with a monogenic causative variant absent") # NO POINTS CAN EXPORT.


##############################################################################################################################
