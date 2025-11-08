#box
library(ggplot2)
library(reshape2)

###Senes_vs_Prolif
#Senes_vs_Prolif
data <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/SP_em_scaled_sig.csv", sep = "\t", row.names = 1)
data$Gene <- rownames(data)
##
#upregulated genes
SP_genes_of_up<- c("C3", "ADAMTS19", "ITIH5", "PDPN", "FAM238C")

# # Filter data
filtered_data <- data[data$Gene %in% SP_genes_of_up, ]

# Convert data from wide format to long format
long_data <- melt(filtered_data, id.vars = "Gene", variable.name = "Condition", value.name = "Expression")

# Draw box plot
ggplot(long_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of Gene Expression for Senes_vs_Prolif", x = "Gene_of_up", y = "Expression")
##
#downregulated genes
SP_genes_of_down<- c("NEK2", "GABRA4", "WNT7B", "DEPDC1-AS1", "LCT-AS1")

# # Filter data
filtered_data <- data[data$Gene %in% SP_genes_of_down, ]

# Convert data from wide format to long format
long_data <- melt(filtered_data, id.vars = "Gene", variable.name = "Condition", value.name = "Expression")

# Draw box plot
ggplot(long_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of Gene Expression for Senes_vs_Prolif", x = "Gene_of_down", y = "Expression")

#Senes_vs_Prolif_adjust
data <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/BSP_em_scaled_sig.csv", sep = "\t", row.names = 1)

data$Gene <- rownames(data)
##
#upregulated genes
ad_SP_genes_of_up<- c("C3", "ADAMTS19", "ITIH5", "PDPN", "FAM238C")

# # Filter data
filtered_data <- data[data$Gene %in% ad_SP_genes_of_up, ]

# Convert data from wide format to long format
long_data <- melt(filtered_data, id.vars = "Gene", variable.name = "Condition", value.name = "Expression")

# Draw box plot
ggplot(long_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of Gene Expression for Senes_vs_Prolif_adjust", x = "Gene_of_up", y = "Expression")
##
#downregulated genes
ad_SP_genes_of_down<- c("NEK2", "GABRA4", "WNT7B", "DEPDC1-AS1", "LCT-AS1")

# # Filter data
filtered_data <- data[data$Gene %in% ad_SP_genes_of_down, ]

# Convert data from wide format to long format
long_data <- melt(filtered_data, id.vars = "Gene", variable.name = "Condition", value.name = "Expression")

# Draw box plot
ggplot(long_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of Gene Expression for Senes_vs_Prolif_adjust", x = "Gene_of_down", y = "Expression")

###
###Senes_MtD_vs_Senes
#Senes_MtD_vs_Senes
data <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/SMS_em_scaled_sig.csv", sep = "\t", row.names = 1)
data$Gene <- rownames(data)
##
#upregulated genes
SMS_genes_of_up<- c("INHBB", "EGLN3", "CXCR4", "OR7E19P", "ENSG00000233818")

# # Filter data
filtered_data <- data[data$Gene %in% SMS_genes_of_up, ]

# Convert data from wide format to long format
long_data <- melt(filtered_data, id.vars = "Gene", variable.name = "Condition", value.name = "Expression")

# Draw box plot
ggplot(long_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of Gene Expression for Senes_MtD_vs_Senes", x = "Gene_of_up", y = "Expression")
##
#downregulated genes
SMS_genes_of_down<- c("MTCO3P12", "MTATP6P2", "CXCL10", "CXCL6", "VCAM1")

# # Filter data
filtered_data <- data[data$Gene %in% SMS_genes_of_down, ]

# Convert data from wide format to long format
long_data <- melt(filtered_data, id.vars = "Gene", variable.name = "Condition", value.name = "Expression")

# Draw box plot
ggplot(long_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of Gene Expression for Senes_MtD_vs_Senes", x = "Gene_of_down", y = "Expression")

#Senes_MtD_vs_Senes_adjust
data <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/BSP_em_scaled_sig.csv", sep = "\t", row.names = 1)

data$Gene <- rownames(data)
##
#upregulated genes
ad_SMS_genes_of_up<- c("INHBB", "EGLN3", "CXCR4", "OR7E19P", "ENSG00000233818")

# # Filter data
filtered_data <- data[data$Gene %in% ad_SMS_genes_of_up, ]

# Convert data from wide format to long format
long_data <- melt(filtered_data, id.vars = "Gene", variable.name = "Condition", value.name = "Expression")

# Draw box plot
ggplot(long_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of Gene for Senes_MtD_vs_Senes_adjust", x = "Gene_of_up", y = "Expression")
##
#downregulated genes
ad_SMS_genes_of_down<- c("MTCO3P12", "MTATP6P2", "CXCL10", "CXCL6", "VCAM1")

# # Filter data
filtered_data <- data[data$Gene %in% ad_SMS_genes_of_down, ]

# Convert data from wide format to long format
long_data <- melt(filtered_data, id.vars = "Gene", variable.name = "Condition", value.name = "Expression")

# Draw box plot
ggplot(long_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of Gene for Senes_MtD_vs_Senes_adjust", x = "Gene_of_down", y = "Expression")


###Senes_MtD_vs_Prolif
#Senes_MtD_vs_Prolif
data <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/SMP_em_scaled_sig.csv", sep = "\t", row.names = 1)
data$Gene <- rownames(data)
##
#upregulated genes
SMP_genes_of_up<- c("EGLN3", "LINC01164", "CXCR4", "OR7E19P", "MRGPRX3")

# # Filter data
filtered_data <- data[data$Gene %in% SMP_genes_of_up, ]

# Convert data from wide format to long format
long_data <- melt(filtered_data, id.vars = "Gene", variable.name = "Condition", value.name = "Expression")

# Draw box plot
ggplot(long_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of Gene Expression for Senes_MtD_vs_Prolif", x = "Gene_of_up", y = "Expression")
##
#downregulated genes
SMP_genes_of_down<- c("KRT19", "CARMN", "IL1RL2", "CIDEA", "CPXM1")

# # Filter data
filtered_data <- data[data$Gene %in% SMP_genes_of_down, ]

# Convert data from wide format to long format
long_data <- melt(filtered_data, id.vars = "Gene", variable.name = "Condition", value.name = "Expression")

# Draw box plot
ggplot(long_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of Gene Expression for Senes_MtD_vs_Prolif", x = "Gene_of_down", y = "Expression")

#Senes_MtD_vs_Prolif_adjust
data <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/BSMP_em_scaled_sig.csv", sep = "\t", row.names = 1)

data$Gene <- rownames(data)
##
#upregulated genes
ad_SMP_genes_of_up<- c("EGLN3", "LINC01164", "CXCR4", "OR7E19P", "MRGPRX3")

# # Filter data
filtered_data <- data[data$Gene %in% ad_SMP_genes_of_up, ]

# Convert data from wide format to long format
long_data <- melt(filtered_data, id.vars = "Gene", variable.name = "Condition", value.name = "Expression")

# Draw box plot
ggplot(long_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of Gene for Senes_MtD_vs_Prolif_adjust", x = "Gene_of_up", y = "Expression")
##
#downregulated genes
ad_SMP_genes_of_down<- c("KRT19", "CARMN", "IL1RL2", "CIDEA", "CPXM1")

# # Filter data
filtered_data <- data[data$Gene %in% ad_SMP_genes_of_down, ]

# Convert data from wide format to long format
long_data <- melt(filtered_data, id.vars = "Gene", variable.name = "Condition", value.name = "Expression")

# Draw box plot
ggplot(long_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of Gene for Senes_MtD_vs_Prolif_adjust", x = "Gene_of_down", y = "Expression")



