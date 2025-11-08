#Step1
#read file
DE_Senes_MtD_vs_Prolif = read.table("/Users/xiaran/Desktop/Due 1.26/assignmentDataset and Description-20240114/DE_Senes_MtD_vs_Prolif.csv", header=TRUE, sep="\t")
DE_Senes_MtD_vs_Senes = read.table("/Users/xiaran/Desktop/Due 1.26/assignmentDataset and Description-20240114/DE_Senes_MtD_vs_Senes.csv", header=TRUE, sep="\t")
DE_Senes_vs_Prolif = read.table("/Users/xiaran/Desktop/Due 1.26/assignmentDataset and Description-20240114/DE_Senes_vs_Prolif.csv", header=TRUE, sep="\t")
em = read.table("/Users/xiaran/Desktop/Due 1.26/assignmentDataset and Description-20240114/EM.csv", header=TRUE, sep="\t")
sample_sheet = read.table("/Users/xiaran/Desktop/Due 1.26/assignmentDataset and Description-20240114/sample_sheet.csv", header=TRUE, sep="\t")
GRCh38 = read.table("/Users/xiaran/Desktop/Due 1.26/assignmentDataset and Description-20240114/Human_Background_GRCh38.p13.csv", header=TRUE, sep="\t")
#basic statistics
library(tidyverse)
# Assume the significance threshold is p-value < 0.05
significant_threshold <- 0.05
#Calculate the number of significantly up-regulated and down-regulated genes
count_significant_genes <- function(data) {
  #Delete rows containing NA
  data <- na.omit(data)
  
  #Count the number of significantly up- and down-regulated genes
  upregulated <- sum(data$p < significant_threshold & data$log2fold > 0)
  downregulated <- sum(data$p < significant_threshold & data$log2fold < 0)
  
  return(list("upregulated" = upregulated, "downregulated" = downregulated))
}

# Apply this function to each data set
stats_MtD_vs_Prolif <- count_significant_genes(DE_Senes_MtD_vs_Prolif)
stats_MtD_vs_Senes <- count_significant_genes(DE_Senes_MtD_vs_Senes)
stats_Senes_vs_Prolif <- count_significant_genes(DE_Senes_vs_Prolif)

# Output results
print(stats_MtD_vs_Prolif)
print(stats_MtD_vs_Senes)
print(stats_Senes_vs_Prolif)
#creat master_temp
master_temp = merge(em, GRCh38, by="ID")

master_Senes_MtDvsP = merge(master_temp, DE_Senes_MtD_vs_Prolif, by="ID")
master_Senes_MtDvsS = merge(master_temp, DE_Senes_MtD_vs_Senes, by="ID")
master_Senes_vsP = merge(master_temp, DE_Senes_vs_Prolif, by="ID")

#Remove column
columns_to_remove_MP = c(5,6,7)
master_Senes_MtD_vs_Prolif = master_Senes_MtDvsP[,-columns_to_remove_MP]

columns_to_remove_MS = c(2,3,4)
master_Senes_MtD_vs_Senes = master_Senes_MtDvsS[,-columns_to_remove_MS]

columns_to_remove_P = c(8,9,10)
master_Senes_vs_Prolif = master_Senes_vsP[,-columns_to_remove_P]

#Step2: Generate the first heat map
#Load necessary libraries
library(ggplot2)
library(reshape2)

#Find non-numeric columns in data frame
non_numeric_columns <- sapply(master_temp, function(x) !is.numeric(x))

# Output the name of the non-numeric column
names(non_numeric_columns)[non_numeric_columns]

# Remove non-numeric columns
master_temp_numeric <- master_temp[, !names(master_temp) %in% c("ID", "SYMBOL", "CHROMOSOME", "BIOTYPE")]

# Make sure that all current data frames are numerical types
str(master_temp_numeric)

# Calculate correlation matrix
correlation_matrix <- cor(master_temp_numeric)

# View the correlation matrix if needed
print(correlation_matrix)


# Melt the correlation matrix into long format
melted_correlation_matrix <- melt(correlation_matrix)

#Create the first heatmap
ggplot(data = melted_correlation_matrix, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.title = element_blank()) 

# Create the second heatmap and add numerical markers
ggplot(data = melted_correlation_matrix, aes(Var1, Var2, fill = value)) +
  geom_tile() + 
  geom_text(aes(label = ifelse(abs(value) >= 0.5, sprintf("%.2f", value), "")),
            size = 3, 
            color = "white", 
            vjust = 0.5, hjust = 0.5) + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name="Pearson Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.text.y = element_text(size = 8), 
        axis.title = element_blank(),
        legend.position = "right") + 
  coord_fixed()
###First heat map of only nine genes
library(ggplot2)
library(reshape2)

#Read expression data
expression_matrix <- read.csv("/Users/xiaran/Desktop/Due 1.26/assignmentDataset and Description-20240114/EM.csv", sep = "\t", row.names = 1)

# Calculate correlation matrix
correlation_matrix <- cor(expression_matrix)

# Convert the correlation matrix to long format for use with ggplot
melted_correlation_matrix <- melt(correlation_matrix)

# Draw heat map
ggplot(data = melted_correlation_matrix, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Sample", title = "Correlation Matrix of Samples")

# # Make a second adjusted heatmap
library(ggplot2)
library(reshape2)
library(corrplot)


# Use the corrplot package to draw a heat map with numerical labels (there is something wrong with the layout of this figure)
corrplot(correlation_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45, 
         addCoef.col = "black",  # Add coefficients
         number.cex = 0.7,       # Control the font size of coefficients
         title = "Correlation Matrix of Samples", 
         cl.lim = c(-1, 1))

# Heatmap after adjusting composition (warning will appear)。

corrplot(correlation_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45, 
         addCoef.col = "black",  
         number.cex = 0.6,       # Control the font size of the coefficient and adjust it appropriately to see if it can solve the problem.
         title = "Correlation Matrix of Samples", 
         cl.lim = c(-1, 1),
         mar = c(0,0,1,0))
#Creat PCA
##Create a PCA scatter plot
library(ggplot2)
expression_matrix <- read.csv("/Users/xiaran/Desktop/Due 1.26/assignmentDataset and Description-20240114/EM.csv", sep = "\t", row.names = 1)

#Transpose the expression matrix so that samples are rows and genes are columns
expression_matrix_t <- t(expression_matrix)

# Perform PCA
pca <- prcomp(expression_matrix_t, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x)
#Create sample group data
sample_groups <- data.frame(Sample_Group = c('Prolif', 'Prolif', 'Prolif', 'Senes', 'Senes', 'Senes', 'Senes_MtD', 'Senes_MtD', 'Senes_MtD'))

# Add sample group data to PCA results
pca_df <- cbind(pca_df, sample_groups)

# Draw PCA scatter plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Sample_Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Sample Groups", x = "Principal Component 1", y = "Principal Component 2")
#mean PCA
library(ggplot2)
library(dplyr)

data = em

data_means <- data %>%
  rowwise() %>%
  mutate(Prolif_mean = mean(c(Prolif_1, Prolif_2, Prolif_3)),
         Senes_mean = mean(c(Senes_1, Senes_2, Senes_3)),
         Senes_MtD_mean = mean(c(Senes_MtD_1, Senes_MtD_2, Senes_MtD_3)))

# Extract the average and transpose
means <- data.frame(t(data_means[, c('Prolif_mean', 'Senes_mean', 'Senes_MtD_mean')]))
# Perform PCA
pca <- prcomp(means, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x)
#Add conditional label
pca_df$Condition <- rownames(pca_df)
# Plot PCA results
ggplot(pca_df, aes(x=PC1, y=PC2, label=Condition, color=Condition)) +
  geom_point(size=4) +
  geom_text(vjust=2, hjust=2) +
  theme_minimal() +
  ggtitle("PCA of Mean Gene Expression Data") +
  xlab("Principal Component 1") +
  ylab("Principal Component 2")

