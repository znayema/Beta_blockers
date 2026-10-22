#script from Christine Lary - edited by Jennifer Spillane on Sept 25, 2023
rm(list=ls())

setwd("/work/larylab/Jennifer/beta_blockers/")
#renv::install("ggplot2")
#renv::install("tidyverse")
#renv::install("dplyr")
library(ggplot2)
library(tidyverse)
library(dplyr)

dxa = read.csv(file="/work/larylab/Mayo/data/HDAC_dxa_fit5.csv",stringsAsFactors = F)
lab = read.csv(file="/work/larylab/Mayo/data/HDAC_lab_fit5.csv",stringsAsFactors = F)
bp = read.csv(file="/work/larylab/Mayo/data/HDAC_bp_fit5.csv",stringsAsFactors = F)
cat = read.csv(file="/work/larylab/Mayo/data/HDAC_cat_fit1.csv",stringsAsFactors = F)
idxa = read.csv(file="/work/larylab/Mayo/data/imp_dxa_fit5.csv",stringsAsFactors = F)
ilab = read.csv(file="/work/larylab/Mayo/data/imp_lab_fit5.csv",stringsAsFactors = F)
ibp = read.csv(file="/work/larylab/Mayo/data/imp_bp_fit5.csv",stringsAsFactors = F)
icat = read.csv(file="/work/larylab/Mayo/data/imp_cat_fit1.csv",stringsAsFactors = F)
hdac4snp = "chr2.239301385.C.G_C" # interacts with spine
adrb1_1 = "chr10.114043616.T.C_T" # interacts with fnbmd, drubmd # our top snp
adrb1_2 = "chr10.114044277.A.G_A" # interacts with fnbmd, drubmd

# look up rsIDs for snps
library(biomaRt)

######################
#not working because the operation times out
################################
## Use the default ENSEMBL Variation Mart & Human dataset
snp = useEnsembl(biomart = "snp", 
                 dataset = "hsapiens_snp")

## Create an example set of coordinates as a dataframe
SNP_M <- data.frame(CHR = c(10,10), START = c(114046473, 114046473), END = c(114046473, 114046473))

## Combine these into the format chr:start:end
## It's important to include the end even if it's a single base, 
## otherwise it searches to the end of the chromosome
coords <- apply(SNP_M, 1, paste, collapse = ":")
coords

######################
#not working because it needs "snp" from above
################################
## Submit the query
getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele'),
      filters = c('chromosomal_region'), 
      values = coords, 
      mart = snp)



dxakeep = dxa %>%
  filter(SNP_pval<0.05 | SNP_AvsPl_pval<0.05 | SNP_NvsPl_pval<0.05 | 
          SNP_P40vsPl_pval<0.05 | SNP_P80vsPl_pval<0.05)
idxakeep = idxa %>%  
  filter(SNP_pval<0.05 | SNP_AvsPl_pval<0.05 | SNP_NvsPl_pval<0.05 | 
          SNP_P40vsPl_pval<0.05 | SNP_P80vsPl_pval<0.05)
######################
#not working because it needs "snp" from above
################################
dxakeep[,c("snp","dxa_var")]
idxakeep[,c("snp","dxa_var")]

makeplot <- function(thesnp,thetrait,df) {
  filename = paste("pharmacogenomics/results",thesnp,"_",thetrait,".tiff",sep="")
  #print(filename)
  cols = c("Pl0_pred","Pl1_pred","Pl2_pred","A0_pred","A1_pred","A2_pred",
           "N0_pred","N1_pred","N2_pred","P40.0_pred","P40.1_pred","P40.2_pred",
           "P80.0_pred","P80.1_pred","P80.2_pred")
  get = df %>%
    filter(snp==thesnp) %>%
    filter(dxa_var==thetrait)
  get_mean = get[,"Mean_BMD0"]
  get = get[,cols]/get_mean
  themeans = data.frame(t(get))
  names(themeans) = "emmean"
  themeans$alleles = c(0,1,2,0,1,2,0,1,2,0,1,2,0,1,2)
  themeans$Group = c("Pl","Pl","Pl","A","A","A","N","N","N","P40","P40","P40","P80","P80","P80")
  themeans$Group = factor(themeans$Group,ordered=T,levels=c("Pl","A","N","P40","P80"))
  tiff(filename)
  ggplot(data=themeans, aes(x=alleles, y=emmean, fill=Group)) +
    geom_bar(stat="identity", position = "dodge") +
    theme_minimal() +
    scale_fill_manual(values=c("lightgrey","steelblue","steelblue1","lightblue2","lightblue3")) +
    theme(text = element_text(size = 16)) +
    xlab("Genotype") + ylab(paste("% Change ",toupper(thetrait)))
  dev.off()
}
df = idxa
thetrait = "fnbmd"
thesnp = "chr10.114043616.T.C_T"
makeplot(thesnp,thetrait,df)
