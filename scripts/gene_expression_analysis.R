## Code for gene expression analysis of ecDNA-associated tumors

library(DESeq2)
library(GSVA)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

########################################################
####### Load data
########################################################
figureDir <- file.path('figures')
dataDir <- file.path('data')
rnaDir <- file.path(dataDir, 'mrna-expression')
vsd <- readRDS(file.path(rnaDir, "Metastatic_Tumors_VST_data.rds"))
dds <- readRDS(file.path(rnaDir, "Metastatic_Tumors_VST_data_ecDNA_Annotation.rds"))

amplicon_df <- read_csv("data/AC_Results_v2.csv")
genes <- scan("data/ecDNA_Genes.txt", what=character())

# Yes=ecDNA, No=Absent
plot_colors <- c("No" = "#BEAED4", "Yes" = "#7FC97F")

########################################################
####### Expression of ecDNA-associated genes
########################################################
# Get ecDNA amplicon genes
gs <- list(
  'CSCC_0005-M1' = scan(file.path(dataDir, 'CSCC_0005-M1_1_AllGenes.txt'), what=character()),
  'CSCC_0010-M1' = scan(file.path(dataDir, 'CSCC_0010-M1_1_AllGenes.txt'), what=character()),
  # 'CSCC_0011-M1' = scan(file.path(dataDir, 'CSCC_0011-M1_1_AllGenes.txt'), what=character()), # no RNA-Seq samples
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

figure1D <- df %>%
  mutate(geneset = str_replace(geneset, "CSCC_", "")) %>%
  ggplot(aes(geneset, score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = sameSample), pch = 21, size = 3) +
  theme_classic() +
  labs(x = "ecDNA Amplicon", y = "Gene Expression Z-Score", 
       fill = "Sample with\necDNA Amplicon") +
  theme(
    axis.text.x = element_text(angle = 45, vjust=0.5),
    text = element_text(size = 14),
  ) +
  guides(fill = "none") +
  scale_fill_manual(values = plot_colors)
figure1D

## Test whether ecDNA+ samples have higher amplicon expression of their
## amplicon genesets than the other samples
t.test(score ~ sameSample, data=df)


########################################################
####### ecDNA- vs ecDNA+ DEG Analysis
########################################################
## Prep data
stopifnot(all(colnames(dds) == colnames(vsd)))
vsd$ecDNA_Status <- dds$ecDNA_Status
colData(dds)$ecDNA_Status <- factor(dds$ecDNA_Status, levels = c("Absent", "Present"))
design(dds) <- ~ ecDNA_Status
dds <- DESeq(dds)


########################################################
####### PCA
########################################################
plot_colors <- c("Absent" = "#BEAED4", "Present" = "#7FC97F")

pca_df <- plotPCA(vsd, intgroup="ecDNA_Status", returnData=TRUE)

figure1E <- pca_df %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill=ecDNA_Status), pch=21, size=3) +
  theme_classic() +
  scale_fill_manual(values=plot_colors) +
  labs(fill = "ecDNA") +
  theme(text = element_text(size=14))
figure1E

########################################################
####### DEGs
########################################################
results_map <- lfcShrink(dds, coef="ecDNA_Status_Present_vs_Absent")
map_df <- data.frame(results_map) %>% 
  as_tibble() %>% 
  mutate(
    abs_log2FoldChange = abs(log2FoldChange),
    Up_in_ecDNA = ifelse(log2FoldChange > 0, "Yes", "No")
  )
map_df$gene <- rownames(results_map)

# Number of significant DEGs
map_df %>%
  filter(
    !is.na(log2FoldChange), 
    !is.na(padj), 
    abs_log2FoldChange > 0.5, 
    padj < 0.05
  ) %>%
  count()

map_df <- map_df %>%
  arrange(desc(abs_log2FoldChange))
write_csv(map_df, "data/ecDNA_Present_vs_Absent_DEGs_MAP.csv")

## PTGS2
gene <- "PTGS2"
xdf <- plotCounts(dds, intgroup="ecDNA_Status", gene=gene, returnData=T)
figure1F <- xdf %>%
  ggplot(aes(ecDNA_Status, log2(count))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = plot_colors) +
  geom_point(aes(fill = ecDNA_Status), pch=21, size=3) +
  theme_classic() +
  labs(x = "ecDNA", y = "Log2 Counts", title=gene) +
  theme(text = element_text(size=14),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none")
figure1F