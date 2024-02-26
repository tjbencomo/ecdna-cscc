## Description: Check if expression of genes on ecDNA amplicons is higher in
## ecDNA samples vs non-ecDNA samples
## Makes figure 1c

library(DESeq2)
library(GSVA)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

figureDir <- file.path('figures')
dataDir <- file.path('data')
rnaDir <- file.path(dataDir, 'mrna-expression')
vsd <- readRDS(file.path(rnaDir, "Metastatic_Tumors_VST_data.rds"))

amplicon_df <- read_csv("data/AC_Results_v2.csv")
genes <- scan("data/ecDNA_Genes.txt", what=character())


# Get ecDNA amplicon genes
gs <- list(
  'CSCC_0005-M1' = scan(file.path(dataDir, 'CSCC_0005-M1_1_AllGenes.txt'), what=character()),
  'CSCC_0010-M1' = scan(file.path(dataDir, 'CSCC_0010-M1_1_AllGenes.txt'), what=character()),
  'CSCC_0134-M1_1' = scan(file.path(dataDir, 'CSCC_0134-M1_1_AllGenes.txt'), what=character()),
  'CSCC_0134-M1_2' = scan(file.path(dataDir, 'CSCC_0134-M1_4_AllGenes.txt'), what=character())
)

# Score samples
es <- gsva(assay(vsd), gs, verbose = TRUE)
esdf <- data.frame(t(es))



# We want to compare samples across amplicon genesets so we apply a Z-score to
# each amplicon score for direct comparison
df <- esdf %>%
  scale() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sampleID") %>%
  pivot_longer(!sampleID, names_to = "geneset", values_to = "score") %>%
  mutate(geneset = str_replace(geneset, '\\.', '-')) %>%
  mutate(fromSample = str_split(geneset, "_1|_2", simplify = T)[, 1]) %>%
  mutate(sameSample = ifelse(sampleID == fromSample, "Yes", "No"))

# Figure 1c plot
ampliconExpressionPlot <- df %>%
  mutate(geneset = str_replace(geneset, "CSCC_", "")) %>%
  ggplot(aes(geneset, score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = sameSample), pch = 21, size = 2.5) +
  theme_classic() +
  labs(x = "Amplicon Gene Set", y = "Gene Expression Z-Score", 
       fill = "Sample with\necDNA Amplicon") +
  theme(text = element_text(size = 14)) +
  guides(fill = "none")
ampliconExpressionPlot

## Test whether ecDNA+ samples have higher amplicon expression of their
## amplicon genesets than the other samples
t.test(score ~ sameSample, data=df)
wilcox.test(score ~ sameSample, data=df)

df %>%
  ggplot(aes(sameSample, score)) +
  geom_boxplot(aes(fill = sameSample), outlier.shape = NA) +
  geom_point(position = 'jitter') +
  theme_classic() +
  guides(fill = "none") +
  labs(x = "Sample with ecDNA Amplicon", y = "Gene Expression Z-Score", fill = "")

## Save figures
ggsave(
  filename = file.path(figureDir, "Figure1_Amplicon_Expression.svg"),
  plot = ampliconExpressionPlot,
  width = 8,
  height = 8
)
saveRDS(ampliconExpressionPlot, file.path(figureDir, "Figure1c.rds"))
