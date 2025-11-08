#MA
#Senes_MtD_vs_Prolif
library(ggplot2)
library(dplyr)
library(ggrepel) 

data <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/Amaster_Senes_MtD_vs_Prolif_na_omit.csv", sep="\t")  

# Calculate A value and M value
data$A_value <- (log2(rowMeans(data[,c("Prolif_1", "Prolif_2", "Prolif_3")])) + 
                   log2(rowMeans(data[,c("Senes_MtD_1", "Senes_MtD_2", "Senes_MtD_3")])))/2
data$M_value <- data$log2fold

# Identify the five most significantly up- and down-regulated genes
top_genes <- data %>%
  arrange(desc(abs(M_value))) %>%
  slice(1:5)

# Separate top 5 upregulated and downregulated genes
top_up_genes <- data %>%
  filter(M_value > 0) %>%
  arrange(desc(M_value)) %>%
  head(5)  # Select top 5 upregulated genes

top_down_genes <- data %>%
  filter(M_value < 0) %>%
  arrange(M_value) %>%
  head(5)  # Select top 5 downregulated genes

# Combine top up and down genes
top_genes <- rbind(top_up_genes, top_down_genes)

# Draw MA chart and add labels using ggrepel with adjusted parameters
ggplot(data, aes(x=A_value, y=M_value)) + 
  geom_point(aes(color = M_value > 0), alpha=0.5) +
  scale_color_manual(values = c("blue", "pink")) +
  geom_text_repel(
    data=top_genes,
    aes(label=SYMBOL),
    box.padding = 0.35,   
    point.padding = 0.5,  
    segment.color = 'grey50'
  ) +
  theme_minimal() +
  labs(x="Average Expression (A value)", y="Log2 Fold Change (M value)", title="Senes_MtD vs Prolif") +
  geom_hline(yintercept=0, color="red", linetype="dashed") +
  theme(legend.position = "none")  



###
#Senes_MtD_vs_Senes
data <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/Amaster_Senes_MtD_vs_Senes_na_omit.csv", sep="\t")  

# Calculate A value and M value
data$A_value <- (log2(rowMeans(data[,c("Senes_1", "Senes_2", "Senes_3")])) + 
                   log2(rowMeans(data[,c("Senes_MtD_1", "Senes_MtD_2", "Senes_MtD_3")])))/2
data$M_value <- data$log2fold

# Identify the five most significantly up- and down-regulated genes
top_genes <- data %>%
  arrange(desc(abs(M_value))) %>%
  slice(1:5)

# Separate top 5 upregulated and downregulated genes
top_up_genes <- data %>%
  filter(M_value > 0) %>%
  arrange(desc(M_value)) %>%
  head(5)  # Select top 5 upregulated genes

top_down_genes <- data %>%
  filter(M_value < 0) %>%
  arrange(M_value) %>%
  head(5)  # Select top 5 downregulated genes

# Combine top up and down genes
top_genes <- rbind(top_up_genes, top_down_genes)


# Draw MA chart and add labels using ggrepel
ggplot(data, aes(x=A_value, y=M_value)) + 
  geom_point(aes(color = M_value > 0), alpha=0.5) +
  scale_color_manual(values = c("blue", "pink")) +
  geom_text_repel(
    data=top_genes,
    aes(label=SYMBOL),
    box.padding = 0.35,   
    point.padding = 0.5,  
    segment.color = 'grey50'  
  ) +
  theme_minimal() +
  labs(x="Average Expression (A value)", y="Log2 Fold Change (M value)", title="Senes_MtD_vs_Senes") +
  geom_hline(yintercept=0, color="red", linetype="dashed") +
  theme(legend.position = "none")  

###
#Senes_vs_Prolif
data <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/Amaster_Senes_vs_Prolif_na_omit.csv", sep="\t")  

# Calculate A value and M value
data$A_value <- (log2(rowMeans(data[,c("Prolif_1", "Prolif_2", "Prolif_3")])) + 
                   log2(rowMeans(data[,c("Senes_1", "Senes_2", "Senes_3")])))/2
data$M_value <- data$log2fold

# Identify the five most significantly up- and down-regulated genes
top_genes <- data %>%
  arrange(desc(abs(M_value))) %>%
  slice(1:5)

# Separate top 5 upregulated and downregulated genes
top_up_genes <- data %>%
  filter(M_value > 0) %>%
  arrange(desc(M_value)) %>%
  head(5)  # Select top 5 upregulated genes

top_down_genes <- data %>%
  filter(M_value < 0) %>%
  arrange(M_value) %>%
  head(5)  # Select top 5 downregulated genes

# Combine top up and down genes
top_genes <- rbind(top_up_genes, top_down_genes)

# Draw MA chart and add labels using ggrepel with adjusted parameters
ggplot(data, aes(x=A_value, y=M_value)) + 
  geom_point(aes(color = M_value > 0), alpha=0.5) +
  scale_color_manual(values = c("blue", "pink")) +
  geom_text_repel(
    data=top_genes,
    aes(label=SYMBOL),
    box.padding = 0.35,   
    point.padding = 0.5,  
    segment.color = 'grey50'
  ) +
  theme_minimal() +
  labs(x="Average Expression (A value)", y="Log2 Fold Change (M value)", title="Senes_vs_Prolif") +
  geom_hline(yintercept=0, color="red", linetype="dashed") +
  theme(legend.position = "none")  

