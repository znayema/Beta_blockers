
library(Hmisc)
library(gplots)
#library(tidyr)
library(dplyr)
#library(tidyverse)
library(MASS)
library(lme4)

library(tidyr)
library(heatmaply)



library(patchwork)

library(cowplot)
library(ggplot2)
library(ggrepel)
library(calibrate)




osteomirs = read.csv(file="/work/larylab/Norway/osteomiRs.csv",stringsAsFactors = F,strip.white = T)
osteoclastmirs = read.csv(file="/work/larylab/Norway/osteoclast_mirs.csv",stringsAsFactors = F,strip.white = T)
mymirs = c(osteomirs$X,osteoclastmirs$miRNA)
write.csv(mymirs,"mymirslist.csv",row.names=F,quote=F)
remove = c("miR-17-5p","miR-214-3p","miR-29b-3p","miR-31-5p","miR-100-5p",
           "miR-155-5p")
mymirs = mymirs[!(mymirs %in% remove)]

mymirs2 = paste("hsa-",gsub("-[1-9]p","",mymirs),sep="")
common_mirs_cross <- intersect(filtered_mirs  , mymirs2)


setwd("/work/larylab/NAYEMA/Beta_blocker_project/mirna/results")
df1<-read.csv("mirna_cmn.csv", header = T)
head(df1)
df2<-read.csv("table for all common disperation.csv", header=T)
head(df2)

specific_mirs <- df2[df2$miRNA %in% c("hsa-miR-375", "hsa-miR-320a"), ]
specific_mirs <- long[long$miRNA %in% c("hsa-miR-375", "hsa-miR-320a"), ]



long
filtered_mirs <- df2$miRNA



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


df2_mod <- filter(df2,Analysis %in% c("Baseline HR vs LR")) 




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
  labs(x = expression("log"[2] * "FC"), y = expression("-log"[10] * "p-value"))+ggtitle('Baseline HR vs LR')



saveRDS(Endpoint, file = "Baseline HR vs LR.rds")



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



threshold_fdr <- 0.05
df2_filtered <- filter(df1, pVal.Cmn < threshold_fdr)
DF<-read.csv("mayo_long_spread1.csv")


data3<-read.csv("/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/mirna_long_0.05.csv")
data3_ordered <- data3[order(data3[, 2]), ]
write.csv(data3_ordered,"/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/mirna_long_0.05.csv",row.names=F)


getwd()

all = read.csv("/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/mirna_expand_cmn.csv")
long = read.csv("/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/mirna_long.csv")
long.sig = read.csv("/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/mirna_long_0.05.csv")
all
long
values_of_interest <- c("trtAtenolol (high responder):time", "trt2Atenolol:time", "trtAtenolol (low responder):time", "trt2Atenolol", "trtAtenolol (high responder)", "trtAtenolol (low responder)")

long.sig <-long.sig[long.sig[, 2] %in% values_of_interest, ]
miRNAs_of_interest <- long.sig$miRNA


long<-long[long[, 2] %in% values_of_interest, ]
miRNAs_of_interest <- long$miRNA



filtered_all <- all[all$miRNA %in% miRNAs_of_interest, ]

all_mine <- all %>%
  filter(miRNA %in% c("hsa-miR-19a-3p","hsa-miR-186-5p"))
long_mine <- long %>%
  filter(miRNA %in% c("hsa-miR-19a-3p","hsa-miR-186-5p")) 

# data wrangling to make plots
raw_counts = all_mine[,7:30]
raw_label = gsub("_.*","",names(raw_counts))
norm_counts = all_mine[,31:54]
norm_label = gsub("norm_","",names(norm_counts))
norm_label = gsub("_.*","",norm_label)
# raw data frame
row_label = paste(substring(all_mine$Analysis,1,8),all_mine$miRNA,sep="_")
rawdf = data.frame(raw_label,t(raw_counts),row.names=NULL)
names(rawdf) = c("Group",row_label)
names(rawdf) = gsub("-","_",names(rawdf))
# normalized data frame
row_label = paste(substring(all_mine$Analysis,1,8),all_mine$miRNA,sep="_")
normdf = data.frame(norm_label,t(norm_counts),row.names=NULL)
names(normdf) = c("Group",row_label)
names(normdf) = gsub("-","_",names(normdf))
# raw wrangling
# 19a
raw19a = rawdf[,c(1,2,4)]
raw19a$ID = c(1:nrow(raw19a))
raw19a <- raw19a %>%
  tidyr::gather("timepoint", "counts", 2:3) %>%
  mutate(timepoint = substring(timepoint,1,8)) %>%
  mutate(Group = case_when(Group=="Atenolol..high.responder." ~ "Aten HR",
                           Group=="Atenolol..low.responder." ~ "Aten LR",
                           Group=="Placebo" ~ "Placebo")) %>%
  arrange(ID)
raw19a$ID = as.character((raw19a$ID))

# For my miRNAs
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
selected_data <- miRNA_19a_3p_data %>%
  mutate(log_count = log(counts))


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
########################################end boxplots

saveRDS(Endpoint, file = "hsa_miR_136_3p.rds")

data_filtered <- aten_hr_data %>%
  filter(Group == "Aten HR")

# Calculate mean values for each miRNA at Baseline and Endpoint
mean_values <- selected_data %>%
  group_by(miRNA, timepoint, Group) %>%
  summarise(mean_counts = mean(counts))


group_colors <- c("Aten HR" = "blue", "Aten LR" = "forestgreen", "Placebo" = "grey87")

plot<-ggplot(mean_values, aes(x = timepoint, y = miRNA, fill = Group)) +
  geom_tile() +
  labs(title = "Heatmap of mean_counts for each miRNA by Group",
       x = "Timepoint",
       y = "miRNA") +
  scale_fill_manual(values = group_colors) +
  theme_minimal()





ggsave('/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/plot.pdf', plot = plot)




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


# Libraries
library(tidyverse)
install.packages("hrbrthemes")
library(hrbrthemes)
library(viridis)
library(plotly)
library(d3heatmap)

# Load data 
data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/multivariate.csv", header=T, sep=";")
colnames(data) <- gsub("\\.", " ", colnames(data))

# Select a few country
data <- data %>% 
  filter(Country %in% c("France", "Sweden", "Italy", "Spain", "England", "Portugal", "Greece", "Peru", "Chile", "Brazil", "Argentina", "Bolivia", "Venezuela", "Australia", "New Zealand", "Fiji", "China", "India", "Thailand", "Afghanistan", "Bangladesh", "United States of America", "Canada", "Burundi", "Angola", "Kenya", "Togo")) %>%
  arrange(Country) %>%
  mutate(Country = factor(Country, Country))


# Matrix format
mat <- data

rownames(mat) <- mat[,1]
mat <- mat %>% dplyr::select(-Country, -Group, -Continent)
mat <- as.matrix(data)

# Heatmap
#d3heatmap(mat, scale="column", dendrogram = "none", width="800px", height="80Opx", colors = "Blues")

library(heatmaply)
p <- heatmaply(result, 
               dendrogram = "none",
               xlab = "", ylab = "", 
               main = "",
               scale = "column",
               margins = c(60,100,40,20),
               grid_color = "white",
               grid_width = 0.00001,
               titleX = FALSE,
               hide_colorbar = TRUE,
               branches_lwd = 0.1,
               label_names = c("Country", "Feature:", "Value"),
               fontsize_row = 5, fontsize_col = 5,
               labCol = colnames(mat),
               labRow = rownames(mat),
               heatmap_layers = theme(axis.line=element_blank())
)


result_matrix <- as.matrix(result[, -1])

# Create heatmap
heatmaply(result,dendrogram = "none",
          xlab = "", ylab = "", 
          main = "",
          scale = "column",
          margins = c(60,100,40,20),
          grid_color = "white",
          grid_width = 0.00001,
          titleX = FALSE,
          hide_colorbar = TRUE,
          branches_lwd = 0.1,
          label_names = c("Country", "Feature:", "Value"),
          fontsize_row = 5, fontsize_col = 5,
          labCol = colnames(mat),
          labRow = rownames(mat),
          heatmap_layers = theme(axis.line=element_blank()))


heatmaply(as.matrix(result[, -1]),  # Convert to matrix and exclude the first column (miRNA)
        Rowv = NA, Colv = NA,     # Disable row and column clustering
        scale = "row",            # Scale rows (genes) for better visualization
        margins = c(5, 10),        # Adjust margins for row and column names
        main = "Heatmap of miRNA Expression",  # Set the title
        xlab = "Sample Groups",    # Label for x-axis
        ylab = "miRNA",            # Label for y-axis
        col = heat.colors(256)     # Choose a color palette
)



result_matrix <- as.matrix(result[, -1])  # Convert to matrix and exclude the first column (miRNA)

# Set row names
rownames(result_matrix) <- result$miRNA

# Create the heatmap
heatmap_obj<-heatmaply(result_matrix,
        Rowv = NA, Colv = NA,     # Disable row and column clustering
        scale = "row",            # Scale rows (miRNA) for better visualization
        margins = c(5, 10),        # Adjust margins for row and column names
        main = "miRNA Expression at Baseline vs Endpoint",  # Set the title
       # Label for x-axis
        ylab = "miRNA",  
       xlab="Groups",# Label for y-axis
        col = heat.colors(256),   # Choose a color palette
        cexCol = 1,cexRow = 0.8,cex.main = 0.7,las = 2,                      # Rotate x-axis labels by 45 degrees
       xaxt = "n" ,
       grid_color = "black",
       grid_width = 0.00001,    fontsize_row = 10, fontsize_col = 10        # Choose a color palette
)

install.packages("orca")

library(htmlwidgets)
saveWidget(heatmap_obj, file= "/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/heatmap_long.html")
ggsave("/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/heatmap_long.png", plot = heatmap_obj)

plotly::export(p=heatmap_obj, file = "work/larylab/NAYEMA/Beta_blocker_project/mirna/results/heatmap_long.png")
install.packages("webshot")

orca::orca(heatmap_obj, file = "work/larylab/NAYEMA/Beta_blocker_project/mirna/results/heatmap_long.png")

# Save the summary row as a CSV file
write.csv(k, "summary_row.csv", row.names = FALSE)
webshot::install_phantomjs()

library(webshot)

plotly::export(heatmap_obj, file = "work/larylab/NAYEMA/Beta_blocker_project/mirna/results/heatmap_long.png")

library(phantomjs) 
  numeric_data <- k[, sapply(k, is.numeric)]  # Extract only numeric columns

# Set row names to mirna names
rownames(numeric_data) <- k$miRNA

data_melt <- melt(numeric_data) 
ggp <- ggplot(data_melt, aes(X1, X2)) +       
  geom_tile(aes(fill = value))

# Now, use superheat
superheat(numeric_data, scale = TRUE,heat.col.scheme = "red",smooth.heat = TRUE)
miRNA_names <- k$miRNA
library(magrittr)
options(repr.plot.width = 5, repr.plot.height = 10)

library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
heatmap(as.matrix(numeric_data), Colv = NA, Rowv = NA,  col= colorRampPalette(brewer.pal(8, "Blues"))(25), xlab="variable", ylab="car", main="heatmap",cexRow=1.5)
heatmap(as.matrix(numeric_data), Colv = NA, Rowv = NA, col = color_palette[color_indices],
        xlab = "variable", ylab = "car", main = "heatmap", cexRow = 1.5)
# Create the heatmap with superheat
heatmap(as.matrix(numeric_data)), scale = "column", col = "blue", margins = c(2, 15), Rowv = NA,
        add.expr = grid(nx = NULL, ny = NULL, col = "black", lty = 1))
data(mtcars)

library(superheat)
superheat(mtcars, scale = TRUE)


# Use grayscale color palette based on differences from mean
color_palette <- color_palette <- colorRampPalette(c("lavender", "purple"))(256)

# Calculate color indices based on differences
color_indices <- findInterval(mean_values, seq(min(mean_values), max(mean_values), length.out = 256))

# Set the size of the plotting device
options(repr.plot.width = 5, repr.plot.height = 10)

# Create a heatmap with grayscale color palette
heatmap(as.matrix(numeric_data), Colv = NA, Rowv = NA, main = "heatmap", cexRow = 1.5,cex.axis = 0.7, cexCol = 1)




library(ggplot2)

# Create a ggplot object
heatmap_plot <- ggplot(numeric_data2, aes(x = factor(c("Baseline", "Endpoint")), y = MiRNA)) +
  
  # Create a filled rectangle (heatmap)
  geom_tile(aes(fill = Endpoint - Baseline), color = "white") +
  
  # Add gradient color scale
  scale_fill_gradient(low = "white", high = "brown") +
  
  # Set axis labels and plot title
  labs(x = "Timepoint", y = "miRNA", title = "Heatmap of Baseline and Endpoint Differences") +
  
  # Customize appearance
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))  # Adjust text size for y-axis

# Display the plot
print(heatmap_plot)











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



# read in list of mirs related to bone
osteomirs = read.csv(file="/work/larylab/Norway/osteomiRs.csv",stringsAsFactors = F,strip.white = T)
osteoclastmirs = read.csv(file="/work/larylab/Norway/osteoclast_mirs.csv",stringsAsFactors = F,strip.white = T)
mymirs = c(osteomirs$X,osteoclastmirs$miRNA)
write.csv(mymirs,"mymirslist.csv",row.names=F,quote=F)
remove = c("miR-17-5p","miR-214-3p","miR-29b-3p","miR-31-5p","miR-100-5p",
           "miR-155-5p") # these duplicate miR without the 3p or 5p
mymirs = mymirs[!(mymirs %in% remove)]
mymirs2 = paste("hsa-",gsub("-[1-9]p","",mymirs),sep="")

# read in results from Mayo study
all = read.csv("/work/larylab/NAYEMA/Beta_blocker_project/mirna/results/mirna_expand_cmn.csv",stringsAsFactors = F)

# next annotate with mymirs
all$mymirs = ""

search = gsub("-[1-9]p","",all$miRNA)
for (i in 1:length(mymirs)) {
  get = which(search %in% mymirs2[i])
  if (length(get)>0){
    all$mymirs[get] = mymirs[i]
  }
}

all_top <- all %>%
  filter(adj.pVal.Cmn<0.5)
write.csv(all_top,"top_all_expand.csv",row.names=F,quote=F)

# total mirs found
# calculate enrichment p-value
all_mirs = unique(c(all$miRNA))
sig_mirs = unique(c(all_top$miRNA))
all_mirs_mine = unique(all$mymirs)
all_mirs_mine = all_mirs_mine[all_mirs_mine != ""]
sig_mirs_mine = unique(c(all_top$mymirs))
sig_mirs_mine = sig_mirs_mine[sig_mirs_mine != ""]

c_all = length(unique(all_mirs_mine))
c_sig = length(unique(sig_mirs_mine))
s_all = length(unique(all_mirs))
s_sig = length(unique(sig_mirs))
mytable = rbind(c(c_sig,s_sig-c_sig),c(c_all-c_sig,s_all-c_all
                                       -(s_sig-c_sig)))
print(mytable)
#test_chisq = chisq.test(mytable)
test_fisher = fisher.test(mytable)
print(test_fisher)
