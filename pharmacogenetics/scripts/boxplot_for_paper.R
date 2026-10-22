####Boxplot figures to be used in the beta blocker pharmacogenetics paper
####Author: Jennifer Spillane
####Created: 11-15-23
####Last Modified: 12-6-23



setwd("/work/larylab/Jennifer/beta_blockers/")
renv::activate()

#libraries used
library(dplyr)
library(ggplot2)
library(stringi)
library(ggpubr)


#(if you need to save your renv lockfile)
#renv::snapshot()

#(if you need to clear out your environment)
#rm(list = ls())


##We're interested in making boxplots of the raw (individual level) data for certain snps of interest in adrb, adrb2, and hdac4.
##I'll need to merge the phenotype data to the snp data for each of the genes, and then calculate percent change from the baseline to the endpoint.
##Then I can make boxplots for each of the allele frequencies for each allele of interest.




#################################################
###Reading in all the data
#######################################

#the phenotype data - this file contains missing values as "." which feels like about the worst idea ever. So we're reading it in and telling R to count those as missing values so they'll show up as "NA" instead.
pheno <- read.csv('/work/larylab/data/beta_blocker_data/pharmacogenetics/larypilot_pheno.csv', na.strings='.')

#the key that matches up the snps with their locations
mayo_snps <- read.csv('/work/larylab/data/beta_blocker_data/pharmacogenetics/mayo_snps_list.csv')

#the three files of allele frequencies
adrb1_raw <- read.csv('/work/larylab/data/beta_blocker_data/pharmacogenetics/plinkdata/chr10.ADRB1.raw', sep = '\t')
adrb2_raw <- read.csv('/work/larylab/data/beta_blocker_data/pharmacogenetics/plinkdata/chr5.ADRB2.raw', sep = '\t')
hdac4_raw <- read.csv('/work/larylab/data/beta_blocker_data/pharmacogenetics/plinkdata/chr2.HDAC4.raw', sep = '\t')


#each of the raw allele frequency files will have to be merged with the phenotype data, so I'm adjusting the formatting of that file to match them

#change the name of the column in the phenotype file so it matches the raw files
pheno <- pheno %>% 
  rename("FID" = id_plasma)


#################################################
###ADRB1
#######################################

#strip off the numbering from the FID in the raw file (so it matches the phenotype file)
adrb1_raw <- adrb1_raw %>%
  mutate(FID = stri_replace_all_regex(FID, '\\d+_', '')) 

#merge them together based on FID
adrb1_pheno <- merge(adrb1_raw, pheno, by = "FID")

#in the adrb1 dataframe(s), the snvs are only listed by location, not name, 
#so I need to rename the relevant columns. (just doing this manually for the snv we are interested in) 
#also calculating a percent change column
#also rounding the allele frequency column because the imputed data aren't whole numbers and converting that number to a factor
adrb1_pheno <- adrb1_pheno %>% 
  rename("rs12414657" = chr10.114043616.T.C_T, "rs1801253" = chr10.114045297.G.C_G) %>% 
  mutate(percent_change_ultradist = ((drubmd20 - drubmd0) / drubmd0) * 100,
         percent_change_ctx = ((ctx20 - ctx0) / ctx0) * 100,
         percent_change_femneck = ((fnbmd20 - fnbmd0 / fnbmd0) * 100)) %>% 
  #round allele freq and convert to factor
  mutate(allele_freq_rs12414657 = factor(round(rs12414657)),
         allele_freq_rs1801253 = factor(round(rs1801253)))

adrb1_ultradist_rs12414657 <- adrb1_pheno %>% 
  select(FID, percent_change_ultradist, allele_freq_rs12414657, treat) %>% 
  filter(treat %in% c("P80", "Pl")) %>% 
  filter(allele_freq_rs12414657 != 0)

adrb1_femneck_rs12414657 <- adrb1_pheno %>% 
  select(FID, percent_change_femneck, allele_freq_rs12414657, treat) %>% 
  filter(treat %in% c("P40", "Pl"))

adrb1_rs1801253 <- adrb1_pheno %>% 
  select(FID, percent_change_ctx, allele_freq_rs1801253, treat) %>% 
  filter(treat %in% c("P40", "Pl")) %>% 
  filter(allele_freq_rs1801253 != 0)



#################################################
###ADRB2
#######################################

#strip off the numbering from the FID in the raw file (so it matches the phenotype file)
adrb2_raw <- adrb2_raw %>%
  mutate(FID = stri_replace_all_regex(FID, '\\d+_', '')) 

#merge them together based on FID
adrb2_pheno <- merge(adrb2_raw, pheno, by = "FID")

#in the adrb2 dataframe(s), the snvs are only listed by location, not name, 
#so I need to rename the relevant columns. (just doing this manually for the snv we are interested in) 
#also calculating a percent change column
#also rounding the allele frequency column because the imputed data aren't whole numbers and converting that number to a factor
adrb2_pheno <- adrb2_pheno %>% 
  rename("rs1042714" = chr5.148826910.G.C_G, "rs1042711" = chr5.148826785.C.T_C, "rs1042717" = chr5.148827083.G.A_G) %>% 
  mutate(percent_change_body = ((bodybmd20 - bodybmd0) / bodybmd0) * 100,
         percent_change_spine = ((sl1l4bmd20 - sl1l4bmd0) / sl1l4bmd0) * 100) %>% 
  #round allele freq and convert to factor
  mutate(allele_freq_rs1042714 = factor(round(rs1042714)),
         allele_freq_rs1042711 = factor(round(rs1042711)),
         allele_freq_rs1042717 = factor(round(rs1042717)))

adrb2_rs1042714 <- adrb2_pheno %>% 
  select(FID, percent_change_body, allele_freq_rs1042714, treat) %>% 
  filter(treat %in% c("P80", "Pl")) 

adrb2_rs1042711 <- adrb2_pheno %>% 
  select(FID, percent_change_spine, allele_freq_rs1042711, treat) %>% 
  filter(treat %in% c("A", "N", "Pl"))

adrb2_rs1042717 <- adrb2_pheno %>% 
  select(FID, percent_change_spine, allele_freq_rs1042717, treat) %>% 
  filter(treat %in% c("A", "N", "P80", "Pl"))


#################################################
###HDAC4
#######################################

#strip off the numbering from the FID in the raw file (so it matches the phenotype file)
hdac4_raw <- hdac4_raw %>%
  mutate(FID = stri_replace_all_regex(FID, '\\d+_', '')) 

#merge them together based on FID
hdac4_pheno <- merge(hdac4_raw, pheno, by = "FID")

#in the hdac4 dataframe(s), the snvs are only listed by location, not name, so I need to rename the relevant columns
#also calculating a percent change column
#also rounding the allele frequency column because the imputed data aren't whole numbers and converting that number to a factor
hdac4_pheno <- hdac4_pheno %>% 
  rename("rs11124190" = chr2.239301385.C.G_C) %>% 
  mutate(spine_percent_change = ((sl1l4bmd20 - sl1l4bmd0) / sl1l4bmd0) * 100) %>% 
  mutate(ctx_percent_change = ((ctx20 - ctx0) / ctx0) * 100) %>% 
  #round allele freq and convert to factor
  mutate(allele_freq = factor(round(rs11124190)))

hdac4_for_spine_plotting <- hdac4_pheno %>% 
  select(FID, spine_percent_change, allele_freq, treat) %>% 
  filter(treat %in% c("A", "Pl")) %>% 
  filter(allele_freq != 0)

hdac4_for_ctx_plotting <- hdac4_pheno %>% 
  select(FID, ctx_percent_change, allele_freq, treat) %>% 
  filter(treat %in% c("P80", "Pl")) %>% 
  filter(allele_freq != 0)


#################################################
###Making the actual plots
#######################################

#boxplot for rs12414657 ultradistal radius
adrb1_rs12414657_ultradist_boxplot <- ggplot(adrb1_ultradist_rs12414657, aes(x = allele_freq_rs12414657, y = percent_change_ultradist, fill = treat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.18)) +
  ggtitle("ADRB1 - rs12414657 - Ultradistal Radius") +
  ylab("Percent Change") + xlab("Allele Frequency (rs12414657)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20), axis.text = element_text(size = 14),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  scale_fill_manual(name = "Treatment",
                    values = c("forestgreen", "gray87"),
                    labels = c("Propranolol 80", "Placebo"))
  
adrb1_rs12414657_ultradist_boxplot
ggsave('pharmacogenetics/results/boxplots/adrb1_rs12414657_ultradist_prop80.png', plot = adrb1_rs12414657_ultradist_boxplot, device = png, width = 7, height = 5)


#boxplot for rs12414657 femoral neck bmd
adrb1_rs12414657_femneck_boxplot <- ggplot(adrb1_femneck_rs12414657, aes(x = allele_freq_rs12414657, y = percent_change_femneck, fill = treat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.18)) +
  ggtitle("ADRB1 - rs12414657 - Femoral Neck BMD") +
  ylab("Percent Change") + xlab("Allele Frequency (rs12414657)") +
  theme_classic() +
  scale_fill_manual(name = "Treatment",
                    values = c("lightgreen", "gray87"),
                    labels = c("Propranolol 40", "Placebo"))

adrb1_rs12414657_femneck_boxplot
ggsave('pharmacogenetics/results/boxplots/adrb1_rs12414657_femneck_prop40.png', plot = adrb1_rs12414657_femneck_boxplot, device = png, width = 7, height = 5)


#boxplot for rs1801253
adrb1_rs1801253_boxplot <- ggplot(adrb1_rs1801253, aes(x = allele_freq_rs1801253, y = percent_change_ctx, fill = treat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.18)) +
  ggtitle("ADRB1 - rs1801253 - CTX") +
  ylab("Percent Change") + xlab("Allele Frequency (rs1801253)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20), axis.text = element_text(size = 14),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  scale_fill_manual(name = "Treatment",
                    values = c("lightgreen", "gray87"),
                    labels = c("Propranolol 40", "Placebo"))

adrb1_rs1801253_boxplot
ggsave('pharmacogenetics/results/boxplots/adrb1_rs1801253_ctx_prop40.png', plot = adrb1_rs1801253_boxplot, device = png, width = 7, height = 5)


#############################
#boxplots for adrb2

#plot for rs1042714
adrb2_rs1042714_boxplot <- ggplot(adrb2_rs1042714, aes(x = allele_freq_rs1042714, y = percent_change_body, fill = treat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.18)) +
  ggtitle("ADRB2 - rs1042714 - Total Body") +
  ylab("Percent Change") + xlab("Allele Frequency (rs1042714)") +
  theme_classic() +
  scale_fill_manual(name = "Treatment",
                    values = c("forestgreen", "gray87"),
                   labels = c("Propranolol 80", "Placebo"))


adrb2_rs1042714_boxplot
ggsave('pharmacogenetics/results/boxplots/adrb2_rs1042714_totalbody_prop80.png', plot = adrb2_rs1042714_boxplot, device = png, width = 7, height = 5)


#plot for rs1042711
adrb2_rs1042711_boxplot <- ggplot(adrb2_rs1042711, aes(x = allele_freq_rs1042711, y = percent_change_spine, fill = treat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +
  ggtitle("ADRB2 - rs1042711 - Lumbar Spine") +
  ylab("Percent Change") + xlab("Allele Frequency (rs1042711)") +
  theme_classic() +
  scale_fill_manual(name = "Treatment",
                    values = c("skyblue1", "royalblue", "gray87"),
                    labels = c("Atenolol", "Nebivolol", "Placebo"))


adrb2_rs1042711_boxplot
ggsave('pharmacogenetics/results/boxplots/adrb2_rs1042711_lumbarspine_aten_nebiv.png', plot = adrb2_rs1042711_boxplot, device = png, width = 7, height = 5)


#plot for rs1042717
adrb2_rs1042717_boxplot <- ggplot(adrb2_rs1042717, aes(x = allele_freq_rs1042717, y = percent_change_spine, fill = treat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.08)) +
  ggtitle("ADRB2 - rs1042717 - Lumbar Spine") +
  ylab("Percent Change") + xlab("Allele Frequency (rs1042717)") +
  theme_classic() +
  scale_fill_manual(name = "Treatment",
                    values = c("skyblue1", "royalblue", "forestgreen", "gray87"),
                    labels = c("Atenolol", "Nebivolol", "Propranolol 80", "Placebo"))


adrb2_rs1042717_boxplot
ggsave('pharmacogenetics/results/boxplots/adrb2_rs1042717_lumbarspine_aten_nebiv_prop80.png', plot = adrb2_rs1042717_boxplot, device = png, width = 7, height = 5)



#############################
#boxplots for hdac4

#boxplot for hdac4 spine
hdac4_spine_boxplot <- ggplot(hdac4_for_spine_plotting, aes(x = allele_freq, y = spine_percent_change, fill = treat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.18)) +
  ggtitle("HDAC4 - rs11124190 - Spine") +
  ylab("Percent Change") + xlab("Allele Frequency (rs11124190)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20), axis.text = element_text(size = 14),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  scale_fill_manual(name = "Treatment",
                    values = c("skyblue1", "gray87"),
                    labels = c("Atenolol", "Placebo"))


hdac4_spine_boxplot
ggsave('pharmacogenetics/results/boxplots/hdac4_rs11124190_spine_aten.png', plot = hdac4_spine_boxplot, device = png, width = 7, height = 5)


#boxplot for hdac4 ctx
hdac4_ctx_boxplot <- ggplot(hdac4_for_ctx_plotting, aes(x = allele_freq, y = ctx_percent_change, fill = treat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.18)) +
  ggtitle("HDAC4 - rs11124190 - CTX") +
  ylab("Percent Change") + xlab("Allele Frequency (rs11124190)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20), axis.text = element_text(size = 14),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  scale_fill_manual(name = "Treatment",
                    values = c("forestgreen", "gray87"),
                    labels = c("Propranolol 80", "Placebo"))

hdac4_ctx_boxplot
ggsave('pharmacogenetics/results/boxplots/hdac4_rs11124190_ctx_prop80.png', plot = hdac4_ctx_boxplot, device = png, width = 7, height = 5)




#If we ever decide to have all the treatments represented, these are the colors that work well for them and for the dots on top of their boxes:
#scale_fill_manual(name = "Treatment",
#                  values = c("skyblue1", "royalblue", "lightgreen", "forestgreen", "gray87"),
#                  labels = c("Atenolol", "Nebivolol", "Propranolol 40", "Propranolol 80", "Placebo")) +
#scale_color_manual(values = c("skyblue3", "royalblue4", "seagreen3", "darkgreen", "gray63"),
#                   guide = "none")


#################################################
###Combining the plots for publication
#######################################

all_plots <- ggarrange(adrb1_rs12414657_ultradist_boxplot, hdac4_spine_boxplot, adrb1_rs1801253_boxplot, hdac4_ctx_boxplot,
                       labels = c("A", "C", "B", "D"),
                       ncol = 2, nrow = 2, align = "hv", font.label = list(size = 20))

all_plots
ggsave('pharmacogenetics/results/boxplots/all_plots.png', plot = all_plots, device = png, width = 14, height = 10)
