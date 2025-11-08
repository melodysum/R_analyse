#volcano
library(ggplot2)
library(dplyr)
library(ggrepel)

#Senes_MtD_vs_Prolif
data <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/SigAmaster_Senes_MtD_vs_Prolif_na_omit_sig.csv", sep="\t")

data$mlog10p <- -log10(data$p)

# Select significantly up- and down-regulated genes
top_upregulated <- data %>% 
  filter(log2fold > 1, p < 0.05) %>% 
  top_n(5, wt=log2fold)
top_downregulated <- data %>% 
  filter(log2fold < -1, p < 0.05) %>% 
  top_n(5, wt=-log2fold)

labels <- rbind(top_upregulated, top_downregulated)

# Draw a volcano diagram and add labels
volcano_plot <- ggplot(data, aes(x=log2fold, y=mlog10p)) +
  geom_point(alpha=0.5) +
  geom_point(data=subset(data, log2fold > 1 & p < 0.05), 
             colour="red", alpha=0.5) +
  geom_point(data=subset(data, log2fold < -1 & p < 0.05), 
             colour="blue", alpha=0.5) +
  geom_text_repel(data=labels, aes(label=SYMBOL), 
                  box.padding = 0.35, 
                  point.padding = 0.5, 
                  segment.color = 'grey50') +
  ggtitle("Volcano Plot - Senes_MtD_vs_Prolif") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 p-value") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# save Picture
ggsave("volcano_plot_Senes_MtD_vs_Prolif.png", volcano_plot, bg = "white")

##
#Senes_MtD_vs_Senes
data <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/SigAmaster_Senes_MtD_vs_Senes_na_omit.csv", sep="\t")

data$mlog10p <- -log10(data$p)

# Select significantly up- and down-regulated genes
top_upregulated <- data %>% 
  filter(log2fold > 1, p < 0.05) %>% 
  top_n(5, wt=log2fold)
top_downregulated <- data %>% 
  filter(log2fold < -1, p < 0.05) %>% 
  top_n(5, wt=-log2fold)

labels <- rbind(top_upregulated, top_downregulated)

# Draw a volcano diagram and add labels
volcano_plot <- ggplot(data, aes(x=log2fold, y=mlog10p)) +
  geom_point(alpha=0.5) +
  geom_point(data=subset(data, log2fold > 1 & p < 0.05), 
             colour="red", alpha=0.5) +
  geom_point(data=subset(data, log2fold < -1 & p < 0.05), 
             colour="blue", alpha=0.5) +
  geom_text_repel(data=labels, aes(label=SYMBOL), 
                  box.padding = 0.35, 
                  point.padding = 0.5, 
                  segment.color = 'grey50') +
  ggtitle("Volcano Plot - Senes_MtD_vs_Senes") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 p-value") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# save Picture
ggsave("volcano_plot_Senes_MtD_vs_Senes.png", volcano_plot, bg = "white")
###
#Senes_vs_Prolif
data <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/SigAmaster_Senes_vs_Prolif_na_omit.csv", sep="\t")

data$mlog10p <- -log10(data$p)

# Select significantly up- and down-regulated genes
top_upregulated <- data %>% 
  filter(log2fold > 1, p < 0.05) %>% 
  top_n(5, wt=log2fold)
top_downregulated <- data %>% 
  filter(log2fold < -1, p < 0.05) %>% 
  top_n(5, wt=-log2fold)

labels <- rbind(top_upregulated, top_downregulated)

# Draw a volcano diagram and add labels
volcano_plot <- ggplot(data, aes(x=log2fold, y=mlog10p)) +
  geom_point(alpha=0.5) +
  geom_point(data=subset(data, log2fold > 1 & p < 0.05), 
             colour="red", alpha=0.5) +
  geom_point(data=subset(data, log2fold < -1 & p < 0.05), 
             colour="blue", alpha=0.5) +
  geom_text_repel(data=labels, aes(label=SYMBOL), 
                  box.padding = 0.35, 
                  point.padding = 0.5, 
                  segment.color = 'grey50') +
  ggtitle("Volcano Plot - Senes_vs_Prolif") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 p-value") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# save Picture
ggsave("volcano_plot_Senes_vs_Prolif.png", volcano_plot, bg = "white")

