#pathway analysis

library(clusterProfiler)
library(org.Hs.eg.db)
#Senes_MtD_vs_Prolif
data_P_SMP <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/Amaster_Senes_MtD_vs_Prolif_na_omit.csv", sep="\t", header=TRUE)

#Extract gene symbols
genes <- unique(data_P_SMP$SYMBOL)

# Perform enrichment analysis
ego_SMP <- enrichGO(gene         = genes,
                OrgDb        = org.Hs.eg.db,
                keyType      = "SYMBOL",
                ont          = "BP",  # Biological Process (BP), Molecular Function (MF), Cellular Component (CC)
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable     = TRUE)

# View the top ten significantly enriched pathways
head(ego_SMP@result, 10)
##

#Senes_MtD_vs_Senes
data_P_SMS <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/Amaster_Senes_MtD_vs_Senes_na_omit.csv", sep="\t", header=TRUE)

#Extract gene symbols
genes <- unique(data_P_SMS$SYMBOL)

# Perform enrichment analysis
ego_SMS <- enrichGO(gene         = genes,
                OrgDb        = org.Hs.eg.db,
                keyType      = "SYMBOL",
                ont          = "BP",  # Biological Process (BP), Molecular Function (MF), Cellular Component (CC)
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable     = TRUE)

# View the top ten significantly enriched pathways
head(ego_SMS@result, 10)
##
#Senes_vs_Prolif
data_P_SP <- read.csv("/Users/xiaran/Desktop/Due 1.26/Sig_gene/Amaster_Senes_vs_Prolif_na_omit.csv", sep="\t", header=TRUE)

#Extract gene symbols
genes <- unique(data_P_SP$SYMBOL)

# Perform enrichment analysis
ego_SP <- enrichGO(gene         = genes,
                OrgDb        = org.Hs.eg.db,
                keyType      = "SYMBOL",
                ont          = "BP",  # Biological Process (BP), Molecular Function (MF), Cellular Component (CC)
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable     = TRUE)
# View the top ten significantly enriched pathways
head(ego_SP@result, 10)
