####Processing data files for the Mayo beta blocker project - pharmacogenetics
####Author: Jennifer Spillane
####Created: 11-20-23
####Last Modified: 11-29-23

#Griffin did this work originally, but it's hard to tell what's going on in some of the scripts
#I'm going to take the parts that are useful and unify it into a single processing script



#############################################################
####get the workspace set up
#################################################

setwd("/work/larylab/Jennifer/beta_blockers/")
renv::activate()

#renv::install("bioc::rtracklayer")

library(dplyr)
library(stringi)


#(if you need to clear out your environment)
#rm(list = ls())


#############################################################
####read in the original data files (from Beth)
#################################################

#for dxa:
dxa <- read.csv(file = "/work/larylab/data/beta_blocker_data/pharmacogenetics/HDAC_dxa_fit5.csv", stringsAsFactors = F)
idxa <- read.csv(file = "/work/larylab/data/beta_blocker_data/pharmacogenetics/imp_dxa_fit5.csv", stringsAsFactors = F)

#for lab:
lab <- read.csv(file = "/work/larylab/data/beta_blocker_data/pharmacogenetics/HDAC_lab_fit5.csv", stringsAsFactors = F)
ilab <- read.csv(file = "/work/larylab/data/beta_blocker_data/pharmacogenetics/imp_lab_fit5.csv", stringsAsFactors = F)

#the lookup table:
#(we're also altering the way the ch38 locations are written so that they only have the number location 
#and not the chromosome info, which is already in another column)
lookup <- read.csv(file = "/work/larylab/data/beta_blocker_data/pharmacogenetics/mayo_fram_snp_master_list.csv") %>% 
  mutate(ch38 = stri_replace_all_regex(ch38, '[cC]hr\\d+:\\s?(\\d+)', '$1'))



#############################################################
####making the data agree so that they can be merged
#################################################

##Starting with dxa and idxa

#remove "GSA" from the rsids in the dxa file when present
#split the snp column into two bits of info divided by the "_" (since we'll need the effect allele later)
cleaned_dxa <- dxa %>% 
  mutate(snp = stri_replace_all_regex(snp, '(GSA.)?(rs\\d+_[ACTG])', '$2')) %>% 
  mutate(rsid = snp) %>% 
  mutate(rsid = stri_replace_all_regex(rsid, '(rs\\d+)_[ACTG]', '$1')) %>% 
  mutate(snp = stri_replace_all_regex(snp, 'rs\\d+_([ACTG])', '$1')) %>% 
  
  #renaming some of the columns to be more human-readable
  rename(Sample = dxa_var, EA = snp) %>%
  mutate(Sample = recode(Sample, bodybmd = 'Total Body', fnbmd = 'Femoral Neck',
                         fhipbmd = 'Femoral Hip', sl1l4bmd = 'Lumbar Spine', 
                         drubmd = 'Ultradistal Radius', dradbmd = 'Distal Radius'),
         #leaving the non-effect allele as a missing value, since it is not included
         NEA = NA,
         Source = 'dxa') %>%
  
  #merge with the lookup file to get the positions of all the alleles (ch37)
  merge(lookup, by = 'rsid', all = T) %>%
  filter(!is.na(Sample)) %>%
  distinct() %>%
  #get rid of the columns we won't be using
  dplyr::select(-written.allele, -A1, -A2, -rsid.allele, -ch38.reference.allele, -ch38.alternate.allele, -X) %>%
  mutate(Chromosome = as.numeric(gsub('chr', '', Chromosome)))


#make four new columns from the snp column, separating its information
cleaned_idxa <- idxa %>%
  mutate(Chromosome = stri_replace_all_regex(snp, 'chr(\\d+)\\.(\\d+)\\.([ACTG])\\.([ACTG]+_[ACTG])', '$1')) %>% 
  mutate(ch38 = as.numeric(stri_replace_all_regex(snp, '(chr\\d+)\\.(\\d+)\\.([ACTG])\\.([ACTG]+_[ACTG])', '$2'))) %>%
  mutate(EA = stri_replace_all_regex(snp, '(chr\\d+)\\.(\\d+)\\.([ACTG])\\.([ACTG]+_[ACTG])', '$3')) %>%
  mutate(NEA = stri_replace_all_regex(snp, '(chr\\d+)\\.(\\d+)\\.([ACTG])\\.([ACTG]+_[ACTG])', '$4')) %>%
  
  #renaming some of the columns to be more human-readable
  rename(Sample = dxa_var) %>%
  mutate(Sample = recode(Sample, bodybmd = 'Total Body', fnbmd = 'Femoral Neck',
                         fhipbmd = 'Femoral Hip', sl1l4bmd = 'Lumbar Spine', 
                         drubmd = 'Ultradistal Radius', dradbmd = 'Distal Radius')) %>% 
  
  #now we'll merge with the lookup table to get ch37 location and rsid
  merge(lookup, by = c('ch38', 'Chromosome'), all = T) %>%
  filter(!is.na(Sample)) %>%
  distinct() %>%
  #removing the columns we don't need
  dplyr::select(-written.allele, -rsid.allele, -A1, -A2, -ch38.reference.allele, -ch38.alternate.allele, -X, -snp) %>%
  #removing the columns that end in "pred_se"
  dplyr::select(-matches('.*pred_se')) %>% 
  mutate(Source = 'imputed dxa')



##Moving on to lab and ilab - should be a very very similar process

#remove "GSA" from the rsids in the dxa file when present
#split the snp column into two bits of info divided by the "_" (since we'll need the effect allele later)
cleaned_lab <- lab %>% 
  mutate(snp = stri_replace_all_regex(snp, '(GSA.)?(rs\\d+_[ACTG])', '$2')) %>% 
  mutate(rsid = snp) %>% 
  mutate(rsid = stri_replace_all_regex(rsid, '(rs\\d+)_[ACTG]', '$1')) %>% 
  mutate(snp = stri_replace_all_regex(snp, 'rs\\d+_([ACTG])', '$1')) %>% 
  
  #renaming some of the columns to be more human-readable
  rename(Sample = lab_var, EA = snp) %>%
  mutate(#leaving the non-effect allele as a missing value, since it is not included
         NEA = NA,
         Source = 'lab') %>%
  
  #merge with the lookup file to get the positions of all the alleles (ch37)
  merge(lookup, by = 'rsid', all = T) %>%
  filter(!is.na(Sample)) %>%
  distinct() %>%
  #get rid of the columns we won't be using
  dplyr::select(-written.allele, -A1, -A2, -rsid.allele, -ch38.reference.allele, -ch38.alternate.allele, -X, -dxa_var) %>%
  mutate(Chromosome = as.numeric(gsub('chr', '', Chromosome)))


#make four new columns from the snp column, separating its information
cleaned_ilab <- ilab %>%
  mutate(Chromosome = stri_replace_all_regex(snp, 'chr(\\d+)\\.(\\d+)\\.([ACTG])\\.([ACTG]+_[ACTG])', '$1')) %>% 
  mutate(ch38 = as.numeric(stri_replace_all_regex(snp, '(chr\\d+)\\.(\\d+)\\.([ACTG])\\.([ACTG]+_[ACTG])', '$2'))) %>%
  mutate(EA = stri_replace_all_regex(snp, '(chr\\d+)\\.(\\d+)\\.([ACTG])\\.([ACTG]+_[ACTG])', '$3')) %>%
  mutate(NEA = stri_replace_all_regex(snp, '(chr\\d+)\\.(\\d+)\\.([ACTG])\\.([ACTG]+_[ACTG])', '$4')) %>%
  
  #renaming some of the columns to be more human-readable
  rename(Sample = lab_var) %>%
  mutate(Sample = recode(Sample, bodybmd = 'Total Body', fnbmd = 'Femoral Neck',
                         fhipbmd = 'Femoral Hip', sl1l4bmd = 'Lumbar Spine', 
                         drubmd = 'Ultradistal Radius', dradbmd = 'Distal Radius')) %>% 
  
  #now we'll merge with the lookup table to get ch37 location and rsid
  merge(lookup, by = c('ch38', 'Chromosome'), all = T) %>%
  filter(!is.na(Sample)) %>%
  distinct() %>%
  #removing the columns we don't need
  dplyr::select(-written.allele, -rsid.allele, -A1, -A2, -ch38.reference.allele, -ch38.alternate.allele, -X, -snp, -dxa_var) %>%
  #removing the columns that end in "pred_se" 
  dplyr::select(-matches('.*pred_se')) %>% 
  mutate(Source = 'imputed dxa')




#############################################################
####merging the data together and getting rid of duplicates
#################################################

#bind with the non-imputed dataframe first, because those rows are our priority
all_dxa <- rbind(cleaned_dxa, cleaned_idxa) %>% #by itself this is 222 obs of 64 variables
  #this distinct() command will get keep the first unique instance of an rsid and sampling site, and get rid of all supsequent ones
  #(this is why it's important to have the non-imputed data first, so all those unique instances are kept)
  distinct(rsid, Sample, .keep_all = TRUE) #with this line it is 162 obs of 64 variables


#same treatment for the lab data as for the dxa (above)
all_lab <- rbind(cleaned_lab, cleaned_ilab) %>% #by itself this is 288 obs of 64 variables
  distinct(rsid, Sample, .keep_all = TRUE) #with this line it is 208 obs of 64 variables


#if all you want is a cleaned and merged version of the data, you can read out now:
#write.csv(all_dxa, 'pharmacogenetics/data/cleaned_merged_dxa.csv')
#write.csv(all_lab, 'pharmacogenetics/data/cleaned_merged_lab.csv')



#############################################################
####parsing the data by drug and selecting the columns we'll need for plotting
#################################################

#a function to parse the data by drug
parse_data_by_drug <- function(df_in, effect_estimate, standard_error, p_val, Drug_Name, Baseline) {
  df_in %>%
    #filter for p_values significant to the specified drug
    #rename generic p_value column to 'p_value' and round p_value to three decimal places
    mutate(p_value = df_in[[p_val]]) %>% 
    mutate(p_value = round(p_value, 3)) %>%
    #rename columns to be human readable and more representative of what they are
    rename('Baseline' = Baseline, 'Value' = effect_estimate, 'Std.Error' = standard_error, 'p' = p_value) %>%
    
    #add column specifying the significant drug for merge
    mutate(Drug = Drug_Name) %>%
    mutate(Value = Value/Baseline) %>%
    
    #select the columns we'll need going forward
    dplyr::select(N, rsid, ch37, Chromosome, Gene, EA, NEA, Sample, Drug, Value, Std.Error, p, Source)
}


#separate BMD (dxa) data by drug for which the effect was significant
atenolol <- parse_data_by_drug(all_dxa, 'SNP_AvsPl_coef', 'SNP_AvsPl_se', 'SNP_AvsPl_pval', 'Atenolol', 'Mean_BMD0')
nebivolol <- parse_data_by_drug(all_dxa, 'SNP_NvsPl_coef', 'SNP_NvsPl_se', 'SNP_NvsPl_pval', 'Nebivolol', 'Mean_BMD0')
propranolol_40 <- parse_data_by_drug(all_dxa, 'SNP_P40vsPl_coef', 'SNP_P40vsPl_se', 'SNP_P40vsPl_pval', 'Propranolol40', 'Mean_BMD0')
propranolol_80 <- parse_data_by_drug(all_dxa, 'SNP_P80vsPl_coef', 'SNP_P80vsPl_se', 'SNP_P80vsPl_pval', 'Propranolol80', 'Mean_BMD0')

#separate Lab data by drug for which the effect was significant
atenolol_lab <- parse_data_by_drug(all_lab, 'SNP_AvsPl_coef', 'SNP_AvsPl_se', 'SNP_AvsPl_pval', 'Atenolol', 'Mean_Lab0')
nebivolol_lab <- parse_data_by_drug(all_lab, 'SNP_NvsPl_coef', 'SNP_NvsPl_se', 'SNP_NvsPl_pval', 'Nebivolol', 'Mean_Lab0')
propranolol_40_lab <- parse_data_by_drug(all_lab, 'SNP_P40vsPl_coef', 'SNP_P40vsPl_se', 'SNP_P40vsPl_pval', 'Propranolol40', 'Mean_Lab0')
propranolol_80_lab <- parse_data_by_drug(all_lab, 'SNP_P80vsPl_coef', 'SNP_P80vsPl_se', 'SNP_P80vsPl_pval', 'Propranolol80', 'Mean_Lab0')


# Re-merge into large data frames containing only results that were significant and the drug they
# were significant for.
dxa_by_drug <- bind_rows(atenolol, nebivolol, propranolol_40, propranolol_80)
lab_by_drug <- bind_rows(atenolol_lab, nebivolol_lab, propranolol_40_lab, propranolol_80_lab)


#write out to csv files
write.csv(dxa_by_drug, 'pharmacogenetics/data/processedBMD.csv')
write.csv(lab_by_drug, 'pharmacogenetics/data/processedLAB.csv')


