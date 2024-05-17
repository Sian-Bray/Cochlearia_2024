# Load necessary library
library(dplyr)

# Load the two TSV files as dataframes
no_salt_top1057 <- read.delim("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/genes_UK_EU_no-salt_15_BPM_top1057.txt", stringsAsFactors = FALSE)
no_salt_top1perc <- read.delim("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/genes_UK_EU_no-salt_15_BPM_top1perc.txt", stringsAsFactors = FALSE)
original_top1perc <- read.delim("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/genes_UK_dips_Vs_tets_15_BPM_1perc.txt", stringsAsFactors = FALSE)

# Remove duplicates in the 'gene' column from each dataframe
no_salt_top1057_unique <- no_salt_top1057 %>% distinct(gene, .keep_all = TRUE)
no_salt_top1perc_unique <- no_salt_top1perc %>% distinct(gene, .keep_all = TRUE)
original_top1perc_unique <- original_top1perc %>% distinct(gene, .keep_all = TRUE)

# Compare the 'gene' column in dataframe 1 with the 'gene' column in dataframe 2 for overlaps
overlap_1per <- inner_join(no_salt_top1perc_unique, original_top1perc_unique, by = "gene")
overlap_1057 <- inner_join(no_salt_top1057_unique, original_top1perc_unique, by = "gene")

# Write the dataframe containing only the overlaps to a new CSV file
write.csv(overlap_1per, "overlap_genes_1perc.csv", row.names = FALSE)
write.csv(overlap_1057, "overlap_genes_1057.csv", row.names = FALSE)



# ChatGPT: "I have two dataframes in tsv format. Please write some R code to make a histogram of the values in column 'FstH' in dataframe 1.  Onto this in red please mark the locations of the 'FstH' values of all the items in the 'genes' column in data frame 2 that are also present in dataframe 1."
library(ggplot2)

# Read data from TSV files into dataframes
df1 <- read.delim("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/genes_UK_dips_Vs_tets_15_BPM_all.txt", stringsAsFactors = FALSE)
df2 <- read.delim("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/genes_UK_EU_no-salt_15_BPM_top1perc.txt", stringsAsFactors = FALSE)

common_genes <- intersect(df1$gene, df2$gene)
df2_common <- df2[df2$gene %in% common_genes, ]

# Plot histogram of 'FstH' from df1
p <- ggplot(df1, aes(x = FstH)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "lightblue", color = "black") +
  labs(x = "Values in FstH", y = "Density", title = "Histogram of FstH in DataFrame 1")

# Overlay points from df2_common onto the histogram
if (nrow(df2_common) > 0) {
  p <- p + geom_point(data = df2_common, aes(x = FstH, y = 0), color = "red", size = 1)
}

# Display the plot
print(p)


library(ggplot2)

# Read data from TSV files into dataframes
df1 <- read.delim("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/genes_UK_dips_Vs_tets_15_BPM_all.txt", stringsAsFactors = FALSE)
df2 <- read.delim("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/genes_UK_EU_no-salt_15_BPM_top1perc.txt", stringsAsFactors = FALSE)

common_genes <- intersect(df1$gene, df2$gene)
df2_common <- df2[df2$gene %in% common_genes, ]

# Calculate threshold value for top 1% of 'FstH' values in df1
threshold <- 0.318464522

# Plot histogram of 'FstH' from df1
p <- ggplot(df1, aes(x = FstH)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "lightblue", color = "black") +
  labs(x = "Values in FstH", y = "Density", title = "Histogram of FstH in DataFrame 1")

# Overlay points from df2_common onto the histogram
if (nrow(df2_common) > 0) {
  p <- p + geom_point(data = df2_common, aes(x = FstH, y = 0), color = "red", size = 1)
}

# Add vertical line for the threshold value
p <- p + geom_vline(xintercept = threshold, color = "black", linetype = "dashed")

# Display the plot
print(p)





# "I have two dataframes. Please use R to plot the 'FstH' values in dataframe 1 as a histogram. Onto this please plot the values of 'FstH' in dataframe 2 as red dots."
library(ggplot2)

# Assuming df1 and df2 are your dataframes with 'FstH' values
df1 <- read.delim("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/genes_UK_dips_Vs_tets_15_BPM_all.txt", stringsAsFactors = FALSE)
df2 <- read.delim("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/values_of_genes_in_no-salt_list_in_original_full_BPM.txt", stringsAsFactors = FALSE)

# Plot histogram of 'FstH' values from dataframe 1 (df1)
p <- ggplot(df1, aes(x = FstH)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "lightblue", color = "black") +
  labs(x = "Values in FstH", y = "Density", title = "Histogram of FstH in DataFrame 1")

# Add red dots for 'FstH' values from dataframe 2 (df2)
if (!is.null(df2$FstH)) {
  p <- p + geom_point(data = df2, aes(x = FstH, y = 0), color = "red", size = 2)
}

# Display the plot
print(p)





# ChatGPT: "Thank you. Could you also perform a statistical test to see if the values of FstH are significantly elevated for items that overlap at the 'gene' column. FstH values should be from dataframe 1."

# Load necessary libraries
library(ggplot2)

# Load the two TSV files as dataframes
df1 <- read.delim("dataframe1.tsv", stringsAsFactors = FALSE)
df2 <- read.delim("dataframe2.tsv", stringsAsFactors = FALSE)

# Extract the values from column 'FstH' in dataframe 1
values_df1 <- df1$FstH

# Extract the overlapping genes from dataframe 1 and dataframe 2
overlap_genes <- intersect(df1$gene, df2$gene)

# Subset dataframe 1 to include only the overlapping genes
df1_overlap <- df1[df1$gene %in% overlap_genes, ]

# Perform Wilcoxon rank sum test
wilcox_test <- wilcox.test(values_df1 ~ ifelse(df1$gene %in% overlap_genes, "Overlap", "Non-Overlap"))

# Print the result of the Wilcoxon rank sum test
print(wilcox_test)







