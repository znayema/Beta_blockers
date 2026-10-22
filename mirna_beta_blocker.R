library(Hmisc)
library(gplots)
library(edgeR)
library(dplyr)
library(tidyverse)
library(flextable)
library(officer)
library(multiMiR)

data = read.delim(file="/work/larylab/Mayo/data/mature_miRNA_expression.txt",stringsAsFactors = F)
head(data)
key = read.delim(file="/work/larylab/Mayo/data/Copy of Lary Study - Key.txt",stringsAsFactors = F)
head(key)
b = key$X[key$Timepoint=="Baseline"]
b
e = key$X[key$Timepoint=="Endpoint"]
e
t = key$treatment[key$Timepoint=="Baseline"] # this works for all 3 data sets; data contains baseline then endpoint for each sample
t

#Base on time point set baseline dataset and endpoint dataset.

baseline = data %>% 
  dplyr::select(all_of(b))
head(baseline)
endpoint = data %>% 
  dplyr::select(all_of(e))
head(endpoint)
# delta = endpoint - baseline
outdf = NULL
analysis = c("Baseline HR vs LR","Endpoint HR vs LR","Endpoint Atenolol vs Placebo")


#set new names and order data 

newnames <- paste(t,names(baseline), sep="_")
newnames
theorder = order(t)
order(t)


#requires the the dataset to contain only the counts with the row names as the gene ids and the column names as the sample ids:
#Baseline data

thedata1 = baseline
names(thedata1) = newnames
thedata1 = thedata1[,theorder]
head(thedata1)


# 1 edgeR for baseline Atenolol high responder vs low responder

genes = data[,1]
counts1 = thedata1
rownames(counts1) = genes
head(counts1)
dim( counts1 )

#Building the edgeR Object

group = t[theorder]
group
group2 <- ifelse(group=="Placebo","Placebo","Atenolol")
group2

cds1 <- DGEList( counts1 , group = group )
names(cds1)
head(cds1$counts)
cds1$samples

#filter out low count reads since it would be impossible to detect differential expression

cds1 <- cds1[rowSums(1e+06 * cds1$counts/expandAsMatrix(cds1$samples$lib.size, dim(cds1)) > 1) >= 3, ]
#cds2 <- cds[rowSums(cds$counts> 1) >= 3, ]
cds1 <- calcNormFactors( cds1 )
head(cds1)
write.csv(file="baseline.csv",cds1$samples)

#Estimating Dispersions

cds1 <- estimateCommonDisp( cds1 )
names( cds1 )
cds1$common.dispersion


#Tagwise dispersions

cds1 <- estimateTagwiseDisp( cds1 )
names(cds1)
summary(cds1$tagwise.dispersion)


#get common disperationand tigwise disperation

de1.cmn <- exactTest( cds1 , dispersion="common", pair = c( "Atenolol (low responder)", "Atenolol (high responder)" ) )
de1.tgw <- exactTest( cds1 , dispersion="tagwise", pair = c( "Atenolol (low responder)", "Atenolol (high responder)" ) )
names(de1.tgw)
de1.tgw$comparison
head(de1.tgw$table)
head(cds1$counts)

# Top tages for tagwise analysis
options( digits = 3 ) # print only 3 digits
topTags( de1.tgw , n = 20 , sort.by = "p.value" ) # top 20 DE genes

# Back to count matrix for tagwise analysis
cds1$counts[ rownames( topTags( de1.tgw , n = 15 )$table ) , ]
# Sort tagwise results by Fold-Change instead of p-value
resultsByFC1.tgw <- topTags( de1.tgw , n = nrow( de1.tgw$table ) , sort.by = "logFC" )$table
head(resultsByFC1.tgw)

# Store full topTags results table
resultsTbl1.cmn <- topTags( de1.cmn , n = nrow( de1.cmn$table ) )$table
resultsTbl1.tgw <- topTags( de1.tgw , n = nrow( de1.tgw$table ) )$table
head(resultsTbl1.tgw)
head(resultsTbl1.cmn)

table1 <- resultsTbl1.cmn
na <- rownames(table1)
table2 <- table1 %>%
  dplyr::select(-logCPM) %>%
  mutate(Analysis = "Baseline HR vs LR") %>%
  mutate(miRNA = na)
table2
names(table2)[2] <- 'p-value'
table2

#filter table based on FDR 

table2 <-table2[, c("miRNA","Analysis","logFC","p-value","FDR")]

table1common <- table2
table1.01 <- table2  %>%
  filter(FDR < 0.1)
table1.05 <- table2 %>%
  filter(FDR < 0.5)
table1.04 <- table2 %>%
  filter(FDR <0.4)
table1.03 <- table2 %>%
  filter(FDR < 0.3)
table1.02 <- table2 %>%
  filter(FDR < 0.2)


# 2 edgeR for endpoint Atenolol high responder vs low responder
#requires the the dataset to contain only the counts with the row names as the gene ids and the column names as the sample ids:
#Baseline data

thedata2 = endpoint
names(thedata2) = newnames
thedata2 = thedata2[,theorder]
head(thedata2)

genes = data[,1]
counts2 = thedata2
rownames(counts2) = genes
head(counts2)
dim( counts2 )



#Building the edgeR Object

cds2 <- DGEList( counts2 , group = group )
names(cds2)
head(cds2$counts)
cds2$samples


#filter out low count reads since it would be impossible to detect differential expression

cds2 <- cds2[rowSums(1e+06 * cds2$counts/expandAsMatrix(cds2$samples$lib.size, dim(cds2)) > 1) >= 3, ]
#cds2 <- cds[rowSums(cds$counts> 1) >= 3, ]
cds2 <- calcNormFactors( cds2 )
write.csv(file="endpoint.csv",cds2$samples)


#Estimating Dispersions

cds2 <- estimateCommonDisp( cds2 )
names( cds2 )
cds2$common.dispersion


#Tagwise dispersions

cds2 <- estimateTagwiseDisp( cds2 )
names(cds2)
summary(cds2$tagwise.dispersion)


de2.cmn <- exactTest( cds2 , dispersion="common", pair = c( "Atenolol (low responder)", "Atenolol (high responder)" ) )
de2.tgw <- exactTest( cds2 , dispersion="tagwise", pair = c( "Atenolol (low responder)", "Atenolol (high responder)" ) )
names(de1.tgw)
de2.tgw$comparison
head(de2.tgw$table)
head(cds2$counts)

# Top tages for tagwise analysis
options( digits = 3 ) # print only 3 digits
topTags( de2.tgw , n = 20 , sort.by = "p.value" ) # top 20 DE genes

# Back to count matrix for tagwise analysis
cds2$counts[ rownames( topTags( de2.tgw , n = 15 )$table ) , ]
# Sort tagwise results by Fold-Change instead of p-value
resultsByFC2.tgw <- topTags( de2.tgw , n = nrow( de2.tgw$table ) , sort.by = "logFC" )$table
head( resultsByFC2.tgw )

# Store full topTags results table
resultsTbl2.cmn <- topTags( de2.cmn , n = nrow( de2.cmn$table ) )$table
resultsTbl2.tgw <- topTags( de2.tgw , n = nrow( de2.tgw$table ) )$table
head(resultsTbl2.tgw)
table3 <- resultsTbl2.cmn
n <- rownames(table3)
table4 <- table3 %>%
  dplyr::select(-logCPM) %>%
  mutate(Analysis = "Endpoint HR vs LR") %>%
  mutate(miRNA = n)
names(table4)[2] <- 'p-value'
table4 <-table4[, c("miRNA","Analysis","logFC","p-value","FDR")]
head(table4) # common disperation 

table2common <- table4
table2.01 <-table4  %>%
  filter(FDR < 0.1)

table2.05 <- table4 %>%
  filter(FDR < 0.5)

table2.04 <- table4 %>%
  filter(FDR < 0.4)

table2.03 <- table4 %>%
  filter(FDR < 0.3)

table2.02 <- table4 %>%
  filter(FDR < 0.2)


#requires the the dataset to contain only the counts with the row names as the gene ids and the column names as the sample ids:
#Baseline data

# 3 edgeR for baseline Atenolol vs ploceble

thedata3 = endpoint
names(thedata3) = newnames
thedata3 = thedata3[,theorder]
head(thedata3)

genes = data[,1]
counts3 = thedata3
rownames(counts3) = genes
head(counts3)
dim( counts3 )

#Building the edgeR Object


cds3 <- DGEList( counts3 , group = group2 )
names(cds3)
head(cds3$counts)
cds3$samples



#filter out low count reads since it would be impossible to detect differential expression

cds3 <- cds3[rowSums(1e+06 * cds3$counts/expandAsMatrix(cds3$samples$lib.size, dim(cds3)) > 1) >= 3, ]
#cds2 <- cds[rowSums(cds$counts> 1) >= 3, ]
cds3 <- calcNormFactors( cds3 )

#Estimating Dispersions

cds3 <- estimateCommonDisp( cds3 )
names( cds3 )
cds3$common.dispersion

#Tagwise dispersions

cds3 <- estimateTagwiseDisp( cds3 )
names(cds3)
summary(cds3$tagwise.dispersion)


de3.cmn <- exactTest( cds3 , dispersion="common", pair = c( "Placebo" , "Atenolol" ) )
de3.tgw <- exactTest( cds3 , dispersion="tagwise", pair = c( "Placebo" , "Atenolol" ) )


# Top tages for tagwise analysis
options( digits = 3 ) # print only 3 digits
topTags( de3.tgw , n = 20 , sort.by = "p.value" ) # top 20 DE genes

# Back to count matrix for tagwise analysis
cds3$counts[ rownames( topTags( de3.tgw , n = 15 )$table ) , ]
# Sort tagwise results by Fold-Change instead of p-value
resultsByFC3.tgw <- topTags( de3.tgw , n = nrow( de3.tgw$table ) , sort.by = "logFC" )$table
head( resultsByFC3.tgw )


# Store full topTags results table
resultsTbl3.cmn <- topTags( de3.cmn , n = nrow( de3.cmn$table ) )$table
resultsTbl3.tgw <- topTags( de3.tgw , n = nrow( de3.tgw$table ) )$table

table5 <- resultsTbl3.cmn
na <- rownames(table5)
table6 <- table5 %>%
  dplyr::select(-logCPM) %>%
  mutate(Analysis = "Endpoint Atenolol vs Placebo") %>%
  mutate(miRNA = na)
names(table6)[2] <- 'p-value'
table6 <-table6[, c("miRNA","Analysis","logFC","p-value","FDR")]

table3common <- table6
table3.01 <-table6  %>%
  filter(FDR < 0.1)

table3.05 <- table6 %>%
  filter(FDR < 0.5)

table3.04 <- table6 %>%
  filter(FDR < 0.4)

table3.03 <- table6 %>%
  filter(FDR < 0.3)

table3.02 <- table6 %>%
  filter(FDR < 0.2)



#Totally common disperation table

all.common <- rbind(table1common,table2common,table3common)
all.common$logFC <-format(round(all.common$logFC,2),nsmall=2)
all.common$FDR <- format(signif(all.common$FDR,2))
all.common$`p-value`<-format(signif(all.common$`p-value`,2))
write.csv(all.common,"table for all common disperation.csv")

ft <- flextable(all.common)
ft <- theme_vanilla(ft)
ft <- color(ft, part = "footer", color = "#666666")
ft <- colformat_double(x = ft,
                       big.mark=",", digits = 2, na_str = "N/A")
ft <- autofit(ft)
doc = read_docx()
doc = body_add_flextable(doc,ft)
print(doc,"table for all common disperation.docx")


#table for FDR < 0.1

fdr0.1 <- rbind(table1.01,table2.01,table3.01)
fdr0.1$logFC <-format(round(fdr0.1$logFC,2),nsmall=2)
fdr0.1$FDR <- format(signif(fdr0.1$FDR,2))
fdr0.1$`p-value`<-format(signif(fdr0.1$`p-value`,2))
write.csv(fdr0.1,"table for common disperation fdr<0.1.csv")




ft <- flextable(fdr0.1)
ft <- theme_vanilla(ft)
ft <- color(ft, part = "footer", color = "#666666")
ft <- colformat_double(x = ft,
                       big.mark=",", digits = 2, na_str = "N/A")
ft <- autofit(ft)
doc = read_docx()
doc = body_add_flextable(doc,ft)
print(doc,"table for all common disperation fdr<0.1.docx")

#table for FDR < 0.2

fdr0.2 <- rbind(table1.02,table2.02,table3.02)
fdr0.2$logFC <-format(round(fdr0.2$logFC,2),nsmall=2)
fdr0.2$FDR <- format(signif(fdr0.2$FDR,2))
fdr0.2$`p-value`<-format(signif(fdr0.2$`p-value`,2))
write.csv(fdr0.2,"table for common disperation fdr<0.2.csv")
fdr0.2

ft <- flextable(fdr0.2)
ft <- theme_vanilla(ft)
ft <- color(ft, part = "footer", color = "#666666")
ft <- colformat_double(x = ft,
                       big.mark=",", digits = 2, na_str = "N/A")
ft <- autofit(ft)
doc = read_docx()
doc = body_add_flextable(doc,ft)
print(doc,"table for all common disperation fdr<0.2.docx")



