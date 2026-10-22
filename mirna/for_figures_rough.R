#Packages required

library(Hmisc)
library(gplots)
#library(tidyr)
library(dplyr)
#library(tidyverse)
library(MASS)
install.packages("lme4")
library(lme4)
library(tidyr)
library(heatmaply)
install.packages("patchwork")
library(patchwork)
install.packages("cowplot")
library(cowplot)
library(ggplot2)
install.packages("ggrepel")
library(ggrepel)
install.packages("calibrate")
library(calibrate)

setwd("/work/larylab/NAYEMA/Beta_blocker_project/mirna/results")

#Read the data 

df1<-read.csv("mirna_cmn.csv", header = T)
head(df1)
df2<-read.csv("table for all common disperation.csv", header=T)
head(df2)
filtered_df <- df2[df2$p.value < 0.05 & df2$FDR <= 0.05, ]
filtered_mirs <- df2$miRNA


#Check for overlap mirnas in osteomirs and osteoclastmirs
osteomirs = read.csv(file="/work/larylab/Norway/osteomiRs.csv",stringsAsFactors = F,strip.white = T)
osteoclastmirs = read.csv(file="/work/larylab/Norway/osteoclast_mirs.csv",stringsAsFactors = F,strip.white = T)
mymirs = c(osteomirs$X,osteoclastmirs$miRNA)
write.csv(mymirs,"mymirslist.csv",row.names=F,quote=F)
remove = c("miR-17-5p","miR-214-3p","miR-29b-3p","miR-31-5p","miR-100-5p",
           "miR-155-5p")
mymirs = mymirs[!(mymirs %in% remove)]
mymirs2 = paste("hsa-",gsub("-[1-9]p","",mymirs),sep="")

#Check overlap between mymir2 and filtered_mirs
common_mirs_cross <- intersect(filtered_mirs  , mymirs2)


# Two common mirnas observed and check their result in cross sectional and longituidinal study
specific_mirs <- df2[df2$miRNA %in% c("hsa-miR-375", "hsa-miR-320a"), ]
specific_mirs <- long[long$miRNA %in% c("hsa-miR-375", "hsa-miR-320a"), ]


# Volcano plot
#Endpoint
df2_mod <- filter(df2,Analysis %in% c("Endpoint Atenolol vs Placebo")) 
df2_mod$diffexpressed <- "NO"
df2_mod$diffexpressed[df2_mod$logFC > 1 & df2_mod$p.value<0.05] <- "UP"
df2_mod$diffexpressed[df2_mod$logFC < -1 & df2_mod$p.value<0.05] <- "DOWN"
df2_mod$delabel <- ifelse(df2_mod$FDR <= 0.05 & df2_mod$miRNA %in% head(df2_mod[order(df2_mod$FDR), "miRNA"], 30), df2_mod$miRNA, NA)

Endpoint<-ggplot(data = df2_mod, aes(x = logFC, y = -log10(p.value), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_text_repel(point.size = NA, size = 4.2, max.overlaps = Inf, colour = "black",
                  segment.ncp = 3, segment.angle = 20) +
  theme_bw(base_size = 14) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 1) +
  scale_color_manual(values = c("purple", "grey", "orange"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_set(theme_classic(base_size = 12)) +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.title.y = element_text(face = "bold", margin = margin(0, 20, 0, 0), size = rel(1), color = 'black'),
        axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20, 0, 0, 0), size = rel(1), color = 'black'),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = expression("log"[2] * "FC"), y = expression("-log"[10] * "p-value"))+ggtitle('Endpoint Atenolol vs Placebo')

saveRDS(Endpoint, file = "Endpoint Atenolol vs Placebo.rds")


#Baseline
df2_mod <- filter(df2,Analysis %in% c("Baseline HR vs LR")) 
df2_mod$diffexpressed <- "NO"
df2_mod$diffexpressed[df2_mod$logFC > 1 & df2_mod$p.value<0.05] <- "UP"
df2_mod$diffexpressed[df2_mod$logFC < -1 & df2_mod$p.value<0.05] <- "DOWN"
df2_mod$delabel <- ifelse(df2_mod$FDR <= 0.05 & df2_mod$miRNA %in% head(df2_mod[order(df2_mod$FDR), "miRNA"], 30), df2_mod$miRNA, NA)

Baseline<-ggplot(data = df2_mod, aes(x = logFC, y = -log10(p.value), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_text_repel(point.size = NA, size = 4.2, max.overlaps = Inf, colour = "black",
                  segment.ncp = 3, segment.angle = 20) +
  theme_bw(base_size = 14) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 1) +
  scale_color_manual(values = c("purple", "grey", "orange"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_set(theme_classic(base_size = 12)) +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.title.y = element_text(face = "bold", margin = margin(0, 20, 0, 0), size = rel(1), color = 'black'),
        axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20, 0, 0, 0), size = rel(1), color = 'black'),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = expression("log"[2] * "FC"), y = expression("-log"[10] * "p-value"))+ggtitle('Baseline HR vs LR')

saveRDS(Baseline, file = "Baseline HR vs LR.rds")

#Endpoint

df2_mod <- filter(df2,Analysis %in% c("Endpoint HR vs LR")) 
df2_mod$diffexpressed <- "NO"
df2_mod$diffexpressed[df2_mod$logFC > 1 & df2_mod$p.value<0.05] <- "UP"
df2_mod$diffexpressed[df2_mod$logFC < -1 & df2_mod$p.value<0.05] <- "DOWN"
df2_mod$delabel <- ifelse(df2_mod$FDR <= 0.05 & df2_mod$miRNA %in% head(df2_mod[order(df2_mod$FDR), "miRNA"], 30), df2_mod$miRNA, NA)

Endpoint<-ggplot(data = df2_mod, aes(x = logFC, y = -log10(p.value), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_text_repel(point.size = NA, size = 4.2, max.overlaps = Inf, colour = "black",
                  segment.ncp = 3, segment.angle = 20) +
  theme_bw(base_size = 14) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 1) +
  scale_color_manual(values = c("purple", "grey", "orange"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_set(theme_classic(base_size = 12)) +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.title.y = element_text(face = "bold", margin = margin(0, 20, 0, 0), size = rel(1), color = 'black'),
        axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20, 0, 0, 0), size = rel(1), color = 'black'),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = expression("log"[2] * "FC"), y = expression("-log"[10] * "p-value"))+ggtitle('Endpoint HR vs LR')

saveRDS(Endpoint, file = "Endpoint HR vs LR.rds")

#Combine the volcano plots 
volcano1<-readRDS("Baseline HR vs LR.rds")
volcano2<-readRDS("Endpoint HR vs LR.rds")
volcano3<-readRDS("Endpoint Atenolol vs Placebo.rds")

plot_grid(volcano1,volcano2,volcano3, labels = c('A', 'B', 'C'))+ plot_annotation(title = 'Cross‐sectional microRNA analysis')  &
  theme(plot.tag = element_text(face = 'bold'), 
        title= element_text(face = 'bold'))

Crosssectional_microRNA <-plot_grid(volcano1,volcano2,volcano3, labels = c('A', 'B', 'C'))+ plot_annotation(title = 'Cross sectional microRNA analysis')  &
  theme(plot.tag = element_text(face = 'bold'), 
        title= element_text(face = 'bold'))

Crosssectional_microRNA 
ggsave("Crosssectional_microRNA.pdf",height=5,width=4.25,dpi=300,units="in",scale=2)



# Heatmap and box plots


#c1 <- c1 %>%
  tidyr::gather(ncol(c1))  %>%
  mutate(Timepoint = substring(Timepoint, 1, 8)) %>%
  mutate(Group = case_when(
    Group == "Atenolol (high responder)" ~ "Aten HR",
    Group == "Atenolol (low responder)" ~ "Aten LR",
    Group == "Placebo" ~ "Placebo"
  )) %>%
  arrange(ID)


# For my miRNAs make the dataset 
  
raw_counts = filtered_all[,7:30]
raw_label = gsub("_.*","",names(raw_counts))
norm_counts = filtered_all[,31:54]
norm_label = gsub("norm_","",names(norm_counts))
norm_label = gsub("_.*","",norm_label)
# raw data frame
row_label = paste(substring(filtered_all$Analysis,1,8),filtered_all$miRNA,sep="_")
rawdf = data.frame(raw_label,t(raw_counts),row.names=NULL)
names(rawdf) = c("Group",row_label)
names(rawdf) = gsub("-","_",names(rawdf))
# normalized data frame
row_label = paste(substring(filtered_all$Analysis,1,8),filtered_all$miRNA,sep="_")
normdf = data.frame(norm_label,t(norm_counts),row.names=NULL)
names(normdf) = c("Group",row_label)
names(normdf) = gsub("-","_",names(normdf))
# raw wrangling

unique_miRNAs <- unique(sub("^Baseline_|^Endpoint_", "", names(rawdf)[-1]))

# Create a vector with the desired column order
column_order <- c("Group", 
                  sapply(unique_miRNAs, function(mirna) c(paste0("Baseline_", mirna), paste0("Endpoint_", mirna))))

# Flatten the nested vector
column_order <- unlist(column_order)

# Reorder the columns
ordered_rawdf <- rawdf[, column_order]

rawdata=ordered_rawdf
rawdata$ID = c(1:nrow(rawdata))

miRNA_columns <- grep("^Baseline_|^Endpoint_", names(rawdata), value = TRUE)

processed_data_list <- list()

for (miRNA_col in miRNA_columns) {
  processed_data <- rawdata %>%
    select(Group, ID, matches(miRNA_col)) %>%
    tidyr::gather("timepoint", "counts", -c(Group, ID)) %>%
    mutate(timepoint = substring(timepoint, 1, 8)) %>%
    mutate(Group = case_when(Group=="Atenolol..high.responder." ~ "Aten HR",
                             Group=="Atenolol..low.responder." ~ "Aten LR",
                             Group=="Placebo" ~ "Placebo"))%>%
    mutate(miRNA = miRNA_col) %>%
    arrange(ID)
  
  # Store the processed data for each miRNA in the list
  processed_data_list[[miRNA_col]] <- processed_data
}

# Combine the processed data for all miRNAs into a single data frame
final_processed_data <- bind_rows(processed_data_list)


# Filter for Aten HR, Aten LR, and Placebo groups
selected_groups <- c("Aten HR", "Aten LR", "Placebo")

selected_data <- final_processed_data %>%
  filter(Group %in% selected_groups)

selected_data$miRNA <- sub(".*hsa", "hsa", selected_data$miRNA)


#Boxplot-use same code for all the desired mirnas

hsa_miR_375_data <- selected_data %>%
  filter(grepl("hsa_miR_375", miRNA))

adrb1_boxplot <-ggplot(hsa_miR_375_data, aes(x = timepoint, y = counts, fill = Group)) +
  geom_boxplot() +geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.18))+
  labs(title = "Box Plot of Aten HR, Aten LR, and Placebo for hsa_miR_375",
       x = "Timepoint",
       y = "Counts") +
  scale_fill_manual(values = c("Aten HR" = "blue", "Aten LR" = "forestgreen", "Placebo" = "grey87"))+theme_classic()+theme(legend.position ="top")


ggsave('/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/hsa_miR_375.png', plot = adrb1_boxplot,width = 7, height = 5)



#Heat map
# Calculate mean values for each miRNA at Baseline and Endpoint
mean_values <- selected_data %>%
  group_by(miRNA, timepoint, Group) %>%
  summarise(mean_counts = mean(counts))

library(dplyr)
library(tidyr)

# Filter the data first for each condition
filtered_data <- mean_values %>%
  filter(Group %in% c("Aten HR", "Aten LR", "Placebo"))

# Create the desired data frame
result <- filtered_data %>%
  pivot_wider(
    names_from = c("timepoint", "Group"),
    values_from = mean_counts,
    names_sep = "_"
  )

column_order <- c("miRNA", 
                  "Baseline_Aten HR", "Endpoint_Aten HR", 
                  "Baseline_Aten LR", "Endpoint_Aten LR", 
                  "Baseline_Placebo", "Endpoint_Placebo")

# Reorder the columns
result <- result %>%
  select(column_order)
write.csv(result,"/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/heatmap.data.csv",row.names=F)


# Libraries for heatmaps
library(tidyverse)
install.packages("hrbrthemes")
library(hrbrthemes)
library(viridis)
library(plotly)
library(d3heatmap)
library(heatmaply)


result_matrix <- as.matrix(result[, -1])
# Set row names
rownames(result_matrix) <- result$miRNA

# Create the heatmap
heatmap_obj<-heatmap(result_matrix,
        Rowv = NA, Colv = NA,     # Disable row and column clustering
        scale = "row",            # Scale rows (miRNA) for better visualization
        margins = c(5, 10),        # Adjust margins for row and column names
        main = "miRNA Expression at Baseline vs Endpoint",  # Set the title
       # Label for x-axis
        ylab = "miRNA",  
       xlab="Groups",# Label for y-axis
        col = heat.colors(256),   # Choose a color palette
        cexCol = .5,cexRow = 0.5,cex.main = 0.7,las = 2,                      # Rotate x-axis labels by 45 degrees
       xaxt = "n" ,grid="black",
    grid_width = 0.000001,    fontsize_row = 10, fontsize_col = 10        # Choose a color palette
)



library(htmlwidgets)
saveWidget(heatmap_obj, file= "/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/heatmap_long.html")

# Hongyus figures code
# 186
raw186 = rawdf[,c(1,3,5)]
raw186$ID = c(1:nrow(raw186))
raw186 <- raw186 %>%
  tidyr::gather("timepoint", "counts", 2:3) %>%
  mutate(timepoint = substring(timepoint,1,8)) %>%
  mutate(Group = case_when(Group=="Atenolol..high.responder." ~ "Aten HR",
                           Group=="Atenolol..low.responder." ~ "Aten LR",
                           Group=="Placebo" ~ "Placebo")) %>%
  arrange(ID)
raw186$ID = as.character((raw186$ID))
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
norm19a$ID = as.character((raw19a$ID))
# 186
norm186 = normdf[,c(1,3,5)]
norm186$ID = c(1:nrow(norm186))
norm186 <- norm186 %>%
  tidyr::gather("timepoint", "counts", 2:3) %>%
  mutate(timepoint = substring(timepoint,1,8)) %>%
  mutate(Group = case_when(Group=="Atenolol..high.responder." ~ "Aten HR",
                           Group=="Atenolol..low.responder." ~ "Aten LR",
                           Group=="Placebo" ~ "Placebo")) %>%
  arrange(ID)
norm186$ID = as.character((raw186$ID))
# raw plot
# 19a

setwd("/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/")
library(ggplot2)
jpeg("raw19a.jpg",quality=100,res=100)
ggplot(raw19a) +
  aes(x = timepoint, y = counts, col = Group) +
  geom_point() +
  geom_line(aes(group = ID),show.legend = F) +
  labs(title = "miR-19a-3p Raw Counts")
dev.off()
# 186
jpeg("raw186.jpg",quality=100,res=100)
ggplot(raw186) +
  aes(x = timepoint, y = counts, col = Group) +
  geom_point() +
  geom_line(aes(group = ID),show.legend = F) +
  labs(title = "miR-186-5p Raw Counts")
dev.off()

# quick model
mymod = lmer(log(counts)~Group*timepoint+(1|ID),data=norm19a)
summary(mymod)
hist(log(norm19a$counts))
norm19a_fc <- norm19a %>%
  tidyr::spread(key=timepoint,value=counts) %>%
  mutate(fc = log(Endpoint) - log(Baseline)) %>%
  mutate(Group = factor(Group,levels=c("Placebo","Aten LR","Aten HR"))) %>%
  mutate(Group2 = factor(ifelse(Group!="Placebo","Atenolol","Placebo"),levels=c("Placebo","Atenolol")))

jpeg("FC_19a.jpg",quality=100,res=100)
g <- ggplot(data=norm19a_fc, aes(x=Group, y = fc, group = Group))
g + geom_boxplot() + xlab("Group") + ylab("log FC") + theme(text = element_text(size = 20))
dev.off()

jpeg("FC_19a_2groups.jpg",quality=100,res=100)
g <- ggplot(data=norm19a_fc, aes(x=Group2, y = fc, group = Group2))
g + geom_boxplot() + xlab("Group") + ylab("log FC") + theme(text = element_text(size = 20))
dev.off()

summary(lm(fc~Group,data=norm19a_fc))
summary(lm(fc~Group2,data=norm19a_fc))
sd = sd(norm19a_fc$fc)
delta= 0.459/sd
power.t.test(n=100,delta=delta,power=)$power
power.t.test(n=50,delta=,sd = sd,power=0.8)$delta
delta = 0.7229/sd
power.t.test(n=25,delta=delta,power=)$power
mean(norm19a_fc$fc[norm19a_fc$Group != "Placebo"])# 0.6807291
mean(norm19a_fc$fc[norm19a_fc$Group == "Placebo"])# 0.2217455
fc = 0.69
power.t.test(n=50,delta=fc, sd=sd,power=)$power # 87% power

# normalized plot
# 19a
jpeg("norm19a.jpg",quality=100,res=100)
ggplot(norm19a) +
  aes(x = timepoint, y = counts, col = Group) +
  geom_point() +
  geom_line(aes(group = ID),show.legend = F) +
  labs(title = "miR-19a-3p Normalized Counts")
dev.off()

# try to change to grey and blue
jpeg("norm19a_color.jpg",quality=100,res=100)
ggplot(norm19a) +
  aes(x = timepoint, y = counts, col = Group) +
  geom_point() +
  scale_fill_manual(values=c("grey", "blue")) +
  geom_line(aes(group = ID),show.legend = F) +
  labs(title = "miR-19a-3p Normalized Counts")
dev.off()

jpeg("norm19a_2groups_color.jpg",quality=100,res=100)
ggplot(norm19a) +
  aes(x = timepoint, y = counts, col = Treatment) +
  #aes(x = timepoint, y = counts, col = c("grey", "blue") + 
  geom_point() +
  geom_line(aes(group = ID),show.legend = F) +
  labs(title = "miR-19a-3p") +
  theme(text = element_text(size = 20)) + 
  scale_color_manual(values = c("#1b98e0","#999999"))
dev.off()
# 186
jpeg("norm186.jpg",quality=100,res=100)
ggplot(norm186) +
  aes(x = timepoint, y = counts, col = Group) +
  geom_point() +
  geom_line(aes(group = ID),show.legend = F) +
  labs(title = "miR-186-5p Normalized Counts")
dev.off()






# Pathway analysis


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multiMiR")





library(multiMiR)

df<-read.csv("mirna_long_0.05.csv")

miRNA_list <- unique(df$miRNA)




db.tables = multimir_dbTables()
db.tables

curr_vers  <- vers_table[1, "VERSION"]


miRNA_list <- c(
  "hsa-miR-1262", "hsa-miR-144-5p", "hsa-miR-4454", "hsa-miR-4999-5p",
  "hsa-miR-6882-3p", "hsa-miR-92a-1-5p", "hsa-miR-136-3p", "hsa-miR-20a-3p",
  "hsa-miR-26b-3p", "hsa-miR-323a-3p", "hsa-miR-4446-3p", "hsa-miR-449a",
  "hsa-miR-485-5p", "hsa-miR-539-3p", "hsa-miR-6882-5p", "hsa-miR-194-3p",
  "hsa-miR-3199", "hsa-miR-873-5p", "hsa-miR-106b-5p", "hsa-miR-1285-5p",
  "hsa-miR-132-3p", "hsa-miR-421", "hsa-miR-452-5p", "hsa-miR-4753-5p",
  "hsa-miR-590-3p", "hsa-miR-1538", "hsa-miR-429", "hsa-miR-874-5p",
  "hsa-miR-2355-3p", "hsa-miR-580-3p", "hsa-miR-6842-5p", "hsa-miR-100-5p",
  "hsa-miR-106a-5p", "hsa-miR-107", "hsa-miR-15b-5p", "hsa-miR-18a-3p",
  "hsa-miR-19a-3p", "hsa-miR-210-3p", "hsa-miR-3187-3p", "hsa-miR-338-3p",
  "hsa-miR-345-5p", "hsa-miR-374a-5p", "hsa-miR-378c", "hsa-miR-425-5p",
  "hsa-miR-4467", "hsa-miR-4654", "hsa-miR-4667-5p", "hsa-miR-542-3p",
  "hsa-miR-548ay-5p", "hsa-miR-548d-5p", "hsa-miR-548n", "hsa-miR-6859-5p",
  "hsa-miR-6862-5p", "hsa-miR-6894-5p", "hsa-miR-885-5p"
)


library(multiMiR)


example <- get_multimir(mirna = miRNA_list, summary = TRUE)

mirna_target_list <- example[, c("mature_mirna_id", "target_symbol")]


write.csv(example, file = "example_data.csv", row.names = FALSE)



mirna_target_data <- data.frame(
  database = example$database,
  mature_mirna_acc = example$mature_mirna_acc,
  mature_mirna_id = example$mature_mirna_id,
  target_symbol = example$target_symbol,
  target_entrez = example$target_entrez,
  target_ensembl = example$target_ensembl,
  experiment = example$experiment,
  support_type = example$support_type,
  pubmed_id = example$pubmed_id
)

# Save the data frame to a CSV file
write.csv(mirna_target_data, file = "example_data.csv", row.names = FALSE)

saveRDS(example, file = "example_data.rds")
# Detailed information of the validated miRNA-target interaction
head(example$validated)
## target genes that are validated by Luciferase assay
lucif <- select(example1, keytype = "type", keys = "disease.drug",
                columns = columns(example1))
lucif[grep("Luciferase", lucif$experiment), ]
example1@summary[example1@summary[,"target_symbol"] == "KRAS",]
## search 'cisplatin' in disease and drug tables in human
example2 <- get_multimir(disease.drug='cisplatin', table='disease.drug')
nrow(example2@data)
head(example2@data)







