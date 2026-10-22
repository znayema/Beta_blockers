####Figures and Tables to be used on the poster and in the paper for the Mayo beta blocker/pharmacogenetics project
####Author: Jennifer Spillane
####Created: 10-30-23
####Last Modified: 11-29-23

#The alleles here are variants of ADRB1, ADRB2, and HDAC4


setwd("/work/larylab/Jennifer/beta_blockers/")
renv::activate()

library(dplyr)
library(flextable)

renv::snapshot()



#(if you need to clear out your environment)
#rm(list = ls())


#################################################
###Table of all SNPs that had a significant effect on either BMD or CTX *for the paper*
#######################################

#Load in results
mayo_bmd <- read.csv('pharmacogenetics/data/processedBMD.csv')
mayo_lab <- read.csv('pharmacogenetics/data/processedLAB.csv')

#Filter the data so that it is just the SNPs that have a significant effect
mayo_data <- filter(mayo_bmd, p < 0.05) %>%
  rbind(filter(mayo_lab, Sample == 'ctx', p < 0.05))

#Since we want confidence intervals in this version, we have some calculating to do on the effect estimates
#First we'll calculate standard deviation
mayo_extra_nums <- mayo_data %>% 
  mutate(ci_upper = Value + (qnorm(0.975) * Std.Error)) %>% mutate(ci_upper = paste0(round(ci_upper, 3))) %>% 
  mutate(ci_lower = Value - (qnorm(0.975) * Std.Error)) %>% mutate(ci_lower = paste0(round(ci_lower, 3))) %>% 
  mutate(CI = paste0(ci_lower, ", ", ci_upper))

#This is a different version of the "plotting_columns", but for the paper, so the formatting is different
plotting_columns <- mayo_extra_nums %>% 
  mutate(p = round(p, digits = 3)) %>% 
  mutate(Value = paste0(round(Value, 3))) %>% 
  mutate(Value = paste0(Value, " (", CI, "; P = ", p, ")")) %>% 
  select(rsid, Gene, EA, ch37, Drug, Sample, Value) %>% 
  arrange(Gene, ch37) %>% 
  rename('Allele' = rsid, 'Alternate Allele' = EA, 'BMD Site' = Sample, 'Effect Estimate (Lower, Upper; P value)' = Value, 'Position' = ch37)

#Building the actual table
set_flextable_defaults(
  font.size = 6, theme_fun = theme_vanilla,
  padding = 6,
  tabcolsep = 0,
  line_spacing = .8,
  text.align = 'center',
  #background.color = "#EFEFEF"
)

mayo_flextable <- flextable(plotting_columns) %>% 
  set_table_properties(align = 'right', layout = 'autofit')


mayo_flextable

save_as_docx(mayo_flextable, path = "pharmacogenetics/results/new_mayo_flextable_for_paper.docx")



