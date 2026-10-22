####Figures and Tables to be used on the poster and in the paper for the Mayo beta blocker/pharmacogenetics project
####Author: Jennifer Spillane
####Created: 9-28-23
####Last Modified: 11-20-23

#The alleles here are variants of ADRB1, ADRB2, and HDAC4
#For the poster, we are primarily concerned with the alleles ADRB1 rs1801252 G and HDAC4 rs11124190 G and the miRNA miR-19a-3p


setwd("/work/larylab/Jennifer/beta_blockers/")
renv::activate()

library(dplyr)
library(ggplot2)
library(ggsignif)
library(scales)
library(ggforce)
library(ggstats)
library(gt)
library(gtExtras)
library(flextable)
library(patchwork)

renv::snapshot()

#################################################
###The three-panel figure - this will make a bunch of plots, but then we can just pull out the ones we are interested in
###This is *for the poster*, not the paper, but it's the same process at first.
#######################################

mayo_emmeans <- read.csv('pharmacogenetics/data/cleaned_imputed_dxa.csv')
#mayo_snps <- read.csv('/work/larylab/Griffin_Stuff/Griffin_Mayo/results/mayo_snps_list.csv')
#I've moved this on to the new "data" directory, so here is the new (more permanent) location
mayo_snps <- read.csv('/work/larylab/data/beta_blocker_data/pharmacogenetics/mayo_snps_list.csv')


parse_data_by_drug <- function(df_in, effect_estimate, standard_error, p_val, Drug_Name) {
  
  # Calculate bounds of confidence intervals
  df_in %>%
    
    # filter for p_values significant to the specified drug and rename generic p_value column to 'p_value' and round p_value to     three decimal places
    
    mutate(p_value = df_in[[p_val]]) %>% filter(p_value < 0.05) %>% mutate(p_value = round(p_value, 3)) %>%
    
    rename('Baseline' = Mean_BMD0, 'Value' = effect_estimate, 'Std.Error' = standard_error, 'p' = p_value) %>%
    
    #mutate(Value = (Value/Baseline) * 100, Std.Error = (Std.Error/Baseline) * 100) %>%
    
    # Select relevant columns for the specified drug
    #select(SNP, Gene, ch37, ch38, BMD_Site, Value, Std.Error, p, Source, REF, ALT) %>%
    
    # Make sure positions are in numeric terms.
    # Added gsub to remove hidden non-numeric characters from ch38
    mutate(ch38 = sub('.*:', '', ch38)) %>%
    mutate(ch37 = as.numeric(ch37), ch38 = as.numeric(gsub("[^0-9.-]", "", ch38))) %>%
    
    
    # Add column specifying the significant drug for merge
    mutate(Drug = Drug_Name, .before = Sample) 
  
}

#Apply the function to the data for each drug
atenolol <- parse_data_by_drug(mayo_emmeans, 'SNP_AvsPl_coef', 'SNP_AvsPl_se', 'SNP_AvsPl_pval', 'Atenolol')
nebivolol <- parse_data_by_drug(mayo_emmeans, 'SNP_NvsPl_coef', 'SNP_NvsPl_se', 'SNP_NvsPl_pval', 'Nebivolol')
propranolol_40 <- parse_data_by_drug(mayo_emmeans, 'SNP_P40vsPl_coef', 'SNP_P40vsPl_se', 'SNP_P40vsPl_pval', 'Propranolol40')
propranolol_80 <- parse_data_by_drug(mayo_emmeans, 'SNP_P80vsPl_coef', 'SNP_P80vsPl_se', 'SNP_P80vsPl_pval', 'Propranolol80')

#bind the results together
filter_for_significance <- bind_rows(atenolol, nebivolol, propranolol_40, propranolol_80)


prebuild_means <- filter_for_significance%>%
  
  # select out the predicted emmeans for each copy of the alternate allele
  select(-matches('se')) %>%
  select(rsid, Sample, matches('pred'), matches('pred_se'), Source, Gene, EA, NEA)%>%
  
  # Long to wide, put data in terms of SNP 
  tidyr::gather(name, EMmean, -rsid, -Sample, -Source, -Gene, -EA, -NEA) %>%
  arrange(rsid)

prebuild_se <- filter_for_significance%>%
  
  # select out the predicted emmeans for each copy of the alternate allele
  select(rsid, Sample, matches('pred_se'), matches('pred_se'), Source, Gene, EA, NEA)%>%
  
  # Long to wide, put data in terms of SNP 
  tidyr::gather(name, Std.Error, -rsid, -Sample, -Source, -Gene, -EA, -NEA) %>%
  arrange(rsid) %>%
  mutate(name = gsub('*.se','',name))


prebuild <- merge(prebuild_means, prebuild_se) %>%
  # Strip non-numeric characters from the allele
  mutate(Allele = gsub("[^0-9.-]", "", name)) %>%
  
  # Only keep substring after '.' (for P40 and P80) 
  mutate(Allele = gsub('.*\\.','',Allele)) %>%
  
  # Name drug variables so they can be plotted
  mutate(drug_var = case_when(startsWith(name, 'Pl') ~ 'Placebo',
                              startsWith(name, 'A') ~ 'Atenolol',
                              startsWith(name, 'N') ~ 'Nebivolol',
                              startsWith(name, 'P40') ~ 'Propranolol 40',
                              startsWith(name, 'P80') ~ 'Propranolol 80')) %>%
  
  mutate(SNP = case_when(startsWith(rsid, 'rs11124190') ~ 'rs11124190',
                         startsWith(rsid, 'rs1336315111') ~ 'rs1336315111',
                         TRUE ~ rsid)) %>%
  
  # ASSUMES MAYO DATA COUNTS REFERENCE ALLELES, remove this and change
  # labels if Beth says otherwise.
  # Change data to effectively count the number of copies of the alternate
  # allele to better match framingham data.
  
  mutate(Allele = case_when(Allele == '2' ~ 0,
                            Allele == '0' ~ 2,
                            Allele == '1' ~ 1))

prebuild


make_plots <- function(emmeans_data, bmd_site_list) {
  
  for(i in 1:length(bmd_site_list)) {
    print(bmd_site_list[[i]])
    filtered_df <- emmeans_data %>%
      filter(Sample == bmd_site_list[[i]])
    
    
    for(j in 1:length(unique(filtered_df$rsid))) {
      print(unique(filtered_df$rsid)[[j]])
      plot_interaction_effects <- filtered_df %>%
        filter(rsid == unique(filtered_df$rsid)[[j]])
      
      # Make Plot
      p <-  ggplot(data = plot_interaction_effects, aes(fill = drug_var, x = Allele, y = EMmean)) +
        geom_col(position = 'dodge', color = 'black') +
        scale_y_continuous(labels = scales::percent) +
        geom_errorbar(aes(ymin = EMmean - 1.96 * Std.Error,
                          ymax = EMmean + 1.96* Std.Error),
                      width=.2,position=position_dodge(.9),
                      show.legend = F) +
        scale_fill_manual(values = c('Placebo' = 'grey','Atenolol' = 'light blue', 'Nebivolol' = 'blue', 'Propranolol 40' = 'light green', 'Propranolol 80' = 'dark green')) +
        labs(#title = paste(plot_interaction_effects$Gene[[1]], 
                           #paste0(plot_interaction_effects$rsid[[1]]),
                           #plot_interaction_effects$Sample[[1]], '(Imputed)'),
             x="Number of Allele Copies", y="Percent Change From Baseline", fill="Treatment") +
        theme(axis.title = element_text(size=26), axis.text = element_text(size = 20),
              legend.title = element_text(size = 16), legend.text = element_text(size = 14))
      
      
      #print(p) 
      ggsave(ggsave(paste0(unique(filtered_df$rsid)[[j]],'_imp_',plot_interaction_effects$Sample[[1]],'.png'),
                    plot = p,
                    device = png,
                    path = '/work/larylab/Jennifer/beta_blockers/pharmacogenetics/results/plots/'))
      
    }
  }
  
}

make_plots(prebuild, list('Femoral Neck', 'Femoral Hip', 'Ultradistal Radius', 'Total Body', 'Lumbar Spine'))







#################################################
###Table of all SNPs that had a significant effect on either BMD or CTX for the poster
#######################################

#Load in results
mayo_bmd <- read.csv('pharmacogenetics/data/processedBMD.csv')
mayo_lab <- read.csv('pharmacogenetics/data/processedLAB.csv')

#Filter the data so that it is just the SNPs that have a significant effect
mayo_data <- filter(mayo_bmd, p < 0.05) %>%
  rbind(filter(mayo_lab, Sample == 'ctx', p < 0.05))

#Isolating the columns we need for the table and putting the entries in order by SNP name (This was for the poster)
plotting_columns <- mayo_data %>%
  select(rsid, Gene, Drug, Sample, Value, p, ch37) %>%
  mutate(Value = paste0(round(Value, 3))) %>%
  arrange(rsid) %>%
  rename('Allele' = rsid, 'BMD Site' = Sample, 'Effect Estimate' = Value, 'P Value' = p, 'Position' = ch37)

#Building the actual table
set_flextable_defaults(
  font.size = 6, theme_fun = theme_vanilla,
  padding = 1,
  tabcolsep = 0,
  line_spacing = .8,
  text.align = 'center',
  #background.color = "#EFEFEF"
)

mayo_flextable <- flextable(plotting_columns) 

mayo_flextable

save_as_docx(mayo_flextable, path = "pharmacogenetics/results/mayo_flextable.docx")


#Griffin's code for a flextable
adrb1_flextable <- flextable(adrb1) %>%
  add_header_row(values = c('Genotype', 'Mayo Results', 'FHS Results'), colwidths = c(5, 6, 4)) %>%
  align(align = 'center', part = 'all', j = 3:8) %>%
  merge_v(j = 1) %>% merge_v(j = 2) %>% merge_v(j = 3) %>% merge_v(j = 4) %>% merge_v(j = 5) %>%
  merge_v(j = 6) %>% merge_v(j = 7) %>% merge_v(j = 12) %>%
  set_caption(caption = 'ADRB1 Chromosome 10') %>%
  theme_vanilla() %>%
  set_table_properties(align = 'right', layout = 'autofit') %>%
  border_outer() %>%
  border_inner_v()







#################################################
###Hongyu's figure from poster - 19a before and after
#######################################

#this is the figure that was in the abstract of the poster initially from Hongyu's work
#it is the 19a miRNA plot that shows the two timepoints and how the miRNA counts changed over that time

# read in results from Mayo study
all = read.csv("/work/larylab/Mayo/mirna_expand_cmn.csv",stringsAsFactors = F)

# first out results for miR-19a-3p 
all_mine <- all %>%
  filter(miRNA %in% c("hsa-miR-19a-3p","hsa-miR-186-5p"))

# data wrangling to make plots
norm_counts = all_mine[,31:54]
norm_label = gsub("norm_","",names(norm_counts))
norm_label = gsub("_.*","",norm_label)

# normalized data frame
row_label = paste(substring(all_mine$Analysis,1,8),all_mine$miRNA,sep="_")
normdf = data.frame(norm_label,t(norm_counts),row.names=NULL)
names(normdf) = c("Group",row_label)
names(normdf) = gsub("-","_",names(normdf))

# normalized wrangling
# 19a
norm19a = normdf[,c(1,2,4)]
norm19a$ID = c(1:nrow(norm19a))
norm19a <- norm19a %>%
  tidyr::gather("timepoint", "counts", 2:3) %>%
  mutate(timepoint = substring(timepoint,1,8)) %>%
  mutate(Group2 = case_when(Group=="Atenolol..high.responder." ~ "Atenolol",
                            Group=="Atenolol..low.responder." ~ "Atenolol",
                            Group=="Placebo" ~ "Placebo")) %>%
  mutate(Treatment = Group2) %>%
  mutate(Group = case_when(Group=="Atenolol..high.responder." ~ "Aten HR",
                           Group=="Atenolol..low.responder." ~ "Aten LR",
                           Group=="Placebo" ~ "Placebo")) %>%
  arrange(ID)


#Making the plot itself
plot_19a <- ggplot(norm19a) +
  aes(x = timepoint, y = counts, col = Treatment) +
  geom_point() +
  geom_line(linewidth = 1.25, aes(group = ID),show.legend = F) +
  ylab("Read Counts") + xlab("Timepoint") + 
  #labs(title = "miR-19a-3p") +
  theme(axis.title = element_text(size=26), axis.text = element_text(size = 20),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14)) + 
  scale_color_manual(values = c("light blue", "grey"))

#Saving it (in the same way Griffin did so that the fonts and quality match)
ggsave(ggsave('mirna/results/norm19a_2groups_color_test.png'),
              plot = plot_19a,
              device = png)



