## Description: Investigate frequency of cSCC ecDNA genes in NatGen samples
## See get_natgen_genes.py for instructions on creating natgen_ecdna_genes.csv.gz
## 41588_2020_678_MOESM2_ESM.xlsx is supplemental files from Kim 2020 Nature Genetics ecDNA paper
## Makes figures 1a and 1b

library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(patchwork)
library(RColorBrewer)


figureDir <- file.path("figures")

# Load data
# Table with all genes found in ecDNA amplicons and the samples they were found in
gene_df <- read_csv("data/natgen_ecdna_genes.csv.gz")
# Metadata about each sample - lineage + sample type etc
sampleInfo <- read_excel("data/41588_2020_678_MOESM2_ESM.xlsx", sheet=3)

# Count number of tumors in each cancer for frequency calculations
cancerecDNATotals <- sampleInfo %>%
  filter(tumor_or_normal == "tumor", sample_classification == "Circular") %>%
  count(lineage, name = "n_total")

# Combine info
df <- gene_df %>%
  inner_join(sampleInfo)

## See number of times each gene is recurrent in ecDNA samples
countdf <- df %>%
  count(gene)

allGenes <- scan("data/ecDNA_Genes.txt", what=character())
oncoGenes <- scan("data/ecDNA_Oncogenes.txt", what=character())

## Get recurrence counts for genes found in our (cSCC) ecDNA amplicons
countdf %>%
  filter(gene %in% allGenes) %>% 
  arrange(desc(n))

## Get recurrence counts for oncogenes found in our (cSCC) ecDNA amplicons
countdf %>%
  filter(gene %in% oncoGenes) %>% 
  arrange(desc(n))


## See frequency of recurrent genes in head and neck tumors
## This just sees how common recurrent genes are
hnscCountDf <- df %>%
  filter(lineage == "Head and Neck", tumor_or_normal == "tumor") %>%
  count(lineage, gene, name = "n_tumors") %>%
  inner_join(cancerecDNATotals) %>%
  mutate(gene_freq = (n_tumors / n_total) * 100)

skcmCountDf <- df %>%
  filter(lineage == "Skin", tumor_or_normal == "tumor") %>%
  count(lineage, gene, name = "n_tumors") %>%
  inner_join(cancerecDNATotals) %>%
  mutate(gene_freq = (n_tumors / n_total) * 100)

sccCountDf <- df %>%
  filter(
    lineage %in% c("Head and Neck", "Lung Squamous cell"), 
    tumor_or_normal == "tumor"
  ) %>%
  count(lineage, gene, name = "n_tumors") %>%
  inner_join(cancerecDNATotals) %>%
  mutate(gene_freq = (n_tumors / n_total) * 100)

## Look at frequency of our ecDNA amplicon genes from cSCC in HNSC
hnscCountDf %>%
  filter(gene %in% allGenes) %>%
  arrange(desc(gene_freq)) %>%
  View()

## Look at frequency of our ecDNA amplicon genes from cSCC in SKCM
skcmCountDf %>%
  filter(gene %in% allGenes) %>%
  arrange(desc(gene_freq)) %>%
  View()

sccCountDf %>%
  filter(gene %in% allGenes) %>%
  arrange(desc(gene_freq)) %>%
  View()


## Show top 10 genes from our ecDNA amplicons that are most recurrent in the
## the NatGen tumors
top10Genes <- countdf %>%
  filter(gene %in% allGenes) %>%
  slice_max(n, n = 10) %>%
  pull(gene)

allGeneColorCount = df %>%
  filter(gene %in% top10Genes) %>%
  distinct(lineage) %>%
  count() %>%
  pull()
getPalette = colorRampPalette(brewer.pal(12, "Set3"))


df %>%
  filter(gene %in% top10Genes) %>%
  ggplot(aes(fct_infreq(gene))) +
  geom_bar(aes(fill = lineage)) +
  theme_classic() +
  labs(x = "", y = "Number of Tumors", fill = "") +
  scale_fill_manual(values = getPalette(allGeneColorCount)) +
  theme(
    text = element_text(size = 14)
  )

## Repeat but only for oncogenes
top10OncoGenes <- countdf %>%
  filter(gene %in% oncoGenes) %>%
  slice_max(n, n = 10) %>%
  pull(gene)

oncogeneColorCount = df %>%
  filter(gene %in% top10OncoGenes) %>%
  distinct(lineage) %>%
  count() %>%
  pull()

oncogeneNatGenCancerTypePlot <- df %>%
  filter(gene %in% top10OncoGenes) %>%
  ggplot(aes(fct_infreq(gene))) +
  geom_bar(aes(fill = lineage)) +
  theme_classic() +
  labs(x = "", y = "Number of Tumors", fill = "") +
  scale_fill_manual(values = getPalette(oncogeneColorCount)) +
  theme(
    text = element_text(size = 14)
  )
oncogeneNatGenCancerTypePlot

# Figure 1a Plot
oncogeneNatGenCancerTypePlotShort <- df %>%
  filter(gene %in% top10OncoGenes, gene %in% c("MYC", "ARNT", "CALR")) %>%
  mutate(
    lineage = case_when(
      lineage == "Lung Adeno" ~ "Lung Adenocarcinoma",
      lineage == "Skin" ~ "Melanoma",
      lineage == "Uterine Corpus Endometrial" ~ "Endometrial",
      TRUE ~ lineage
    )
  ) %>%
  ggplot(aes(fct_infreq(gene))) +
  geom_bar(aes(fill = lineage)) +
  theme_classic() +
  labs(x = "Oncogene Amplicons", y = "Number of Tumors", fill = "") +
  scale_fill_manual(values = getPalette(oncogeneColorCount)) +
  theme(
    text = element_text(size = 14)
  )
oncogeneNatGenCancerTypePlotShort

## Code to plot probability of detecting ecDNA gene in >1 of our samples
p <- seq(0, .1, .01)
# pNoRecurDetect = pbinom(1, 24, prob=p)
pNoRecurDetect = pbinom(1, 6, prob=p)
pdf <- data.frame(
  prevalence = p,
  probDetectRecurrence = 1 - pNoRecurDetect
)

# Figure 1b plot
probDetectionPlot <- pdf %>%
  ggplot(aes(prevalence * 100, probDetectRecurrence)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(
    x = "Frequency in TCGA/PCAWG Tumors",
    y = "Probability of Detecting Gene in > 1 Tumor"
  ) +
  theme(text = element_text(size = 14)) +
  scale_y_continuous(limits = c(0, .2), breaks = seq(0, .2, by = .05)) +
  scale_x_reverse()
probDetectionPlot

comboPlot <- oncogeneNatGenCancerTypePlotShort | probDetectionPlot

## Save figures
ggsave(
  filename = file.path(figureDir, "Figure1_Oncogenes_NatGen_Frequency.svg"),
  plot = oncogeneNatGenCancerTypePlot,
  width = 8,
  height = 6
)
ggsave(
  filename = file.path(figureDir, "Figure1_Detection_Probability_Plot.svg"),
  plot = probDetectionPlot,
  width = 8, 
  height = 8
)
ggsave(
  filename = file.path(figureDir, "Figure1_Combo.svg"),
  plot = comboPlot,
  width = 16,
  height = 8
)
saveRDS(comboPlot, file.path(figureDir, "Figure1ab.rds"))


