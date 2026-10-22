
install.packages("Hmisc")
library(Hmisc)
library(gplots)
library(tidyr)
library(dplyr)
library(tidyverse)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.12")
#biocLite("edgeR")
library(edgeR)
#library(ptmixed)
library(MASS)
library(lme4)
setwd("/work/larylab/NAYEMA/Beta_blocker_project/mirna_beta_blockers/mirna/results/")

data = read.delim(file="/work/larylab/Mayo/data/mature_miRNA_expression.txt",stringsAsFactors = F)
#themiRNA = data$Mature.miRNA
#themiRNA = substring(themiRNA,nchar(themiRNA)-3,nchar(themiRNA))
#table(themiRNA)
genes = data[,1] # genes names in original data frame "data"
key = read.delim(file="/work/larylab/Mayo/data/Copy of Lary Study - Key.txt",stringsAsFactors = F)
# add additional columns to key
key$ID = sort(rep(1:24,2))
key$time = ifelse(key$Timepoint=="Baseline",1,2)
sampID = names(data)[2:ncol(data)]
sampID = as.numeric(gsub("DM","",sampID))
# order data columns to be the same as key
data = data[,c(1,order(sampID)+1)] # data columns now in same order as key
# now we will order by trt
trt_order = order(key$treatment)
trt_order
key = key[trt_order,]
key
data = data[,c(1,trt_order+1)]
cbind(key$X,names(data)[2:ncol(data)]) # double check key and data still in same order
# set up three data sets: baseline, endpoint, and delta
b = key$X[key$Timepoint=="Baseline"]
e = key$X[key$Timepoint=="Endpoint"]
t = key$treatment[key$Timepoint=="Baseline"] # this works for all 3 data sets; data contains baseline then endpoint for each sample
id = key$ID[key$Timepoint=="Baseline"]
newnames = paste(t,id,sep="_")
baseline = data %>% 
  dplyr::select(all_of(b))
endpoint = data %>% 
  dplyr::select(all_of(e))
#delta = endpoint - baseline
outdf = NULL
analysis = c("Baseline HR vs LR","Endpoint HR vs LR","Endpoint HR vs Placebo","Endpoint LR vs Placebo","Endpoint Atenolol vs Placebo")
for (l in 1:5) {
  if (l == 1) {
    thedata = baseline
  } else {
    thedata = endpoint
  }
  names(thedata) = newnames
  # next load into edgeR
  counts = thedata
  rownames(counts) = genes
  group = t#[theorder]
  group2 = ifelse(group=="Placebo","Placebo","Atenolol")
  if (l<5) {
    cds <- DGEList( counts , group = group )
  } else {
    cds <- DGEList( counts , group = group2 )
  }
  # change from 10 to one for expression filter - back to default - 12/7/2021
  # one read per million in at least 3 samples or 1 RPKM
  cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
  cds <- calcNormFactors( cds )
  if (l==1) write.csv(file="baseline.csv",cds$samples)
  if (l==2) write.csv(file="endpoint.csv",cds$samples)
  cds <- estimateCommonDisp( cds )
  cds <- estimateTagwiseDisp( cds )
  if (l < 3) {
    if (l==1) name = "baseline" else name = "endpoint"
  }
  if (l < 3) {
    de1.cmn <- exactTest( cds , dispersion="common", pair = c( "Atenolol (low responder)", "Atenolol (high responder)" ) )
  } else {
    if (l==3) {
      de1.cmn <- exactTest( cds , dispersion="common", pair = c( "Placebo" , "Atenolol (high responder)" ) )
    } else {
      if (l==4) {
        de1.cmn <- exactTest( cds , dispersion="common", pair = c( "Placebo" , "Atenolol (low responder)" ) )
      } else {
        de1.cmn <- exactTest( cds , dispersion="common", pair = c( "Placebo" , "Atenolol" ) )
      }
    }
  }
  options( digits = 3 ) # print only 3 digits
  # Store full topTags results table
  resultsTbl.cmn <- topTags( de1.cmn , n = nrow( de1.cmn$table ) )$table
  de1.genes.cmn <- rownames( resultsTbl.cmn )[ resultsTbl.cmn$PValue <= 0.05 ]
  colnames( resultsTbl.cmn ) <- c( "logFC" , "logCPM", "pVal.Cmn","adj.pVal.Cmn")
  # Below provides the info to re-order the count matrix to be in line with the order of the results
  wh.rows.cmn <- match( rownames( resultsTbl.cmn ) , rownames( cds$counts ) )
  # Common Results
  combResults.cmn <- cbind( resultsTbl.cmn ,
                            "Cmn.Disp" = cds$common.dispersion ,
                            "UpDown.Cmn" = decideTestsDGE( de1.cmn , p.value = 0.05, adjust.method = "none" )[ wh.rows.cmn,1 ] ,
                            cds$counts[ wh.rows.cmn , ] )
  normcounts = combResults.cmn[,7:ncol(combResults.cmn)]*cds$samples$norm.factors
  if (l==1) {
    normcounts_b = normcounts
  } else {
    if (l==2) {
      normcounts_e = normcounts
    }
  }
  
}

b = normcounts_b
b$miRNA = row.names(b)
e = normcounts_e
e$miRNA = row.names(e)

merged = merge(b,e,by="miRNA",all.x=F,all.y=F)
mirna <- merged$miRNA
mirna

data = read.delim(file="/work/larylab/Mayo/data/mature_miRNA_expression.txt",stringsAsFactors = F)
#themiRNA = data$Mature.miRNA
#themiRNA = substring(themiRNA,nchar(themiRNA)-3,nchar(themiRNA))
#table(themiRNA)
data
genes = data[,1] # genes names in original data frame "data"
key = read.delim(file="/work/larylab/Mayo/data/Copy of Lary Study - Key.txt",stringsAsFactors = F)
key
ex = read.csv(file = "/work/larylab/Mayo/data/larystudy_extravars.csv")
ex
ex1 <- cbind(key,ex)
ex1


data1 <- data %>%
  filter(Mature.miRNA %in% mirna)
data1


data1 <- as.data.frame(t(data1))
data1
colnames(data1) <- data1[1,]

data1 <- data1 %>%
  mutate(X = rownames(data1))

data1 <- data1[2:49,]
data1
arrange(data1, X)

ex1
data1
ex2 <- ex1 %>%
  mutate(time = ifelse(Timepoint == "Baseline",1,2))

ex2
c <- merge(ex2,data1,by="X")
c1 <- c %>%
  mutate(id = gsub("BB","",c$study_id))
c1
as.factor(c1$id)

str(c1)

  c1 
  c1 <- c1 %>%
    mutate(Group = case_when(
      group == "Atenolol (high responder)" ~ "Aten HR",
      group == "Atenolol (low responder)" ~ "Aten LR",
      group == "Placebo" ~ "Placebo"
    )) %>%
    mutate(study_id = as.character(study_id)) %>%
    arrange(study_id) 
  
  
  unique(c1$Timepoint)
  
  baseline_data <- data.frame(
    study_id = unique(ex2$study_id),
    baseline_ctx = ex2[ex2$Timepoint == "Baseline", "ctx"],
    baseline_systolic = ex2[ex2$Timepoint == "Baseline", "systolic"],
    baseline_diastolic = ex2[ex2$Timepoint == "Baseline", "diastolic"],
    baseline_pulse = ex2[ex2$Timepoint == "Baseline", "pulse"]
  )
  
  endpoint_data <- data.frame(
    study_id = unique(ex2$study_id),
    endpoint_ctx = ex2[ex2$Timepoint == "Endpoint", "ctx"],
    endpoint_systolic = ex2[ex2$Timepoint == "Endpoint", "systolic"],
    endpoint_diastolic = ex2[ex2$Timepoint == "Endpoint", "diastolic"],
    endpoint_pulse = ex2[ex2$Timepoint == "Endpoint", "pulse"]
  )
  
  
  
  combined_data$delta_ctx <- combined_data$baseline_ctx - combined_data$endpoint_ctx
  combined_data$delta_systolic <- combined_data$baseline_systolic - combined_data$endpoint_systolic
  combined_data$delta_diastolic <- combined_data$baseline_diastolic - combined_data$endpoint_diastolic
  combined_data$delta_pulse <- combined_data$baseline_pulse - combined_data$endpoint_pulse
  
  correlation_ctx_systolic <- cor.test(combined_data$delta_ctx, combined_data$delta_systolic, method = "pearson")
  correlation_ctx_diastolic <- cor.test(combined_data$delta_ctx, combined_data$delta_diastolic, method = "pearson")
  correlation_ctx_pulse <- cor.test(combined_data$delta_ctx, combined_data$delta_pulse, method = "pearson")
  
  # Create a data frame with results
  correlation_results <- data.frame(
    Variable = c("delta_ctx vs delta_systolic", "delta_ctx vs delta_diastolic", "delta_ctx vs delta_pulse"),
    Correlation = c(correlation_ctx_systolic$estimate, correlation_ctx_diastolic$estimate, correlation_ctx_pulse$estimate),
    P_value = c(correlation_ctx_systolic$p.value, correlation_ctx_diastolic$p.value, correlation_ctx_pulse$p.value),
    Method = c(correlation_ctx_systolic$method, correlation_ctx_diastolic$method, correlation_ctx_pulse$method),
    Conf_Int_Lower = c(correlation_ctx_systolic$conf.int[1], correlation_ctx_diastolic$conf.int[1], correlation_ctx_pulse$conf.int[1]),
    Conf_Int_Upper = c(correlation_ctx_systolic$conf.int[2], correlation_ctx_diastolic$conf.int[2], correlation_ctx_pulse$conf.int[2]),
    Statistic = c(correlation_ctx_systolic$statistic, correlation_ctx_diastolic$statistic, correlation_ctx_pulse$statistic),
    DF = c(correlation_ctx_systolic$parameter, correlation_ctx_diastolic$parameter, correlation_ctx_pulse$parameter)
  )
  
  combined_data
  
  
  write.csv( correlation_results, file="/work/larylab/NAYEMA/Beta_blocker_project/mirna_beta_blockers/mirna/results/correlation_ctx_cardio_overall.csv")
  
  
  library(dplyr)
  
  # Arrange ex2 based on baseline_data structure
  arranged_ex2 <- ex2 %>%
    arrange(study_id, Timepoint) 
  
  ex2$pulse
  
  # Split the arranged_ex2 data into baseline and endpoint data
  baseline_data <- filter(arranged_ex2, Timepoint == "Baseline") %>%
    rename(baseline_ctx = ctx,
           baseline_systolic = systolic,
           baseline_diastolic = diastolic,
           baseline_pulse = pulse) 
  
  endpoint_data <- filter(arranged_ex2, Timepoint == "Endpoint") %>%
    rename(endpoint_ctx = ctx,
           endpoint_systolic = systolic,
           endpoint_diastolic = diastolic,
           endpoint_pulse = pulse) 
  
  combined_data <- cbind(baseline_data, endpoint_data)
  
  combined_data$delta_ctx <- combined_data$baseline_ctx - combined_data$endpoint_ctx
  combined_data$delta_systolic <- combined_data$baseline_systolic - combined_data$endpoint_systolic
  combined_data$delta_diastolic <- combined_data$baseline_diastolic - combined_data$endpoint_diastolic
  combined_data$delta_pulse <- combined_data$baseline_pulse - combined_data$endpoint_pulse
  
  atenolol_low_responder <- combined_data[combined_data$treatment == "Atenolol (low responder)", ]
  atenolol_high_responder <- combined_data[combined_data$treatment == "Atenolol (high responder)", ]
  placebo <- combined_data[combined_data$treatment == "Placebo", ]
  
  # Calculate correlations for each group
  calculate_correlation_tests <- function(data, variable1, variable2, group) {
  cor_test <- cor.test(data[[variable1]], data[[variable2]])
  return(data.frame(
    Variable1 = variable1,
    Variable2 = variable2,
    Correlation = cor_test$estimate,
    P_value = cor_test$p.value,
    Conf_Int_Lower = cor_test$conf.int[1],
    Conf_Int_Upper = cor_test$conf.int[2],
    Test_Statistic = cor_test$statistic,
    Degrees_Freedom = cor_test$parameter,
    Group = group
  ))
}

# Calculate correlation tests for each group and variable pair
cor_tests_atenolol_low <- lapply(c("delta_diastolic", "delta_pulse","delta_systolic"), function(variable) {
  calculate_correlation_tests(atenolol_low_responder, "delta_ctx", variable, "Atenolol (low responder)")
})

cor_tests_atenolol_high <- lapply(c("delta_diastolic", "delta_pulse","delta_systolic"), function(variable) {
  calculate_correlation_tests(atenolol_high_responder, "delta_ctx", variable, "Atenolol (high responder)")
})

cor_tests_placebo <- lapply(c("delta_diastolic", "delta_pulse","delta_systolic"), function(variable) {
  calculate_correlation_tests(placebo, "delta_ctx", variable, "Placebo")
})

# Combine the results into a table
result_table <- do.call(rbind, lapply(c(cor_tests_atenolol_low, cor_tests_atenolol_high, cor_tests_placebo), as.data.frame))

write.csv(combined_data, file="/work/larylab/NAYEMA/Beta_blocker_project/mirna_beta_blockers/mirna/results/data_for_correlation_ctx_cardio.csv")

combined_data$treatment_factor<-as.factor(combined_data$treatment)

a<-lm(delta_ctx~treatment_factor, data = combined_data)

b<-lm(delta_systolic~treatment_factor, data = combined_data)

c<-lm(delta_diastolic~treatment_factor,data=combined_data)
d<-lm(delta_pulse~treatment_factor,data=combined_data)
summary(a)

f <- lm(endpoint_ctx ~ baseline_ctx + treatment_factor, data = combined_data)
summary(f)
