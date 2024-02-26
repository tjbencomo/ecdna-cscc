## Description: Make DESeq and VST objects for further expression analysis
##

library(DESeq2)
library(readr)
library(dplyr)


dataDir <- file.path("data")
rnaDir <- file.path(dataDir, "mrna-expression")

sampleInfo <- read_csv(file.path(dataDir, "Sample_Metadata.csv"))
metTumors <- sampleInfo %>% filter(sample_type == "Metastatic Tumor")

counts_df <- read_tsv(file.path(rnaDir, "unnormalisedGeneCounts.txt"))
colnames(counts_df)[1] <- c("gene")

countsMat <- counts_df %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

metTumorsWithRNA <- metTumors %>%
  filter(SampleName %in% colnames(counts_df))
countsMat <- countsMat[, metTumorsWithRNA$SampleName]

stopifnot(all(metTumorsWithRNA$SampleName == colnames(countsMat)))

dds <- DESeqDataSetFromMatrix(
  countData = countsMat, colData = metTumorsWithRNA, design = ~ 1
)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

saveRDS(dds, file.path(rnaDir, "Metastatic_Tumors_deseq_obj.rds"))
saveRDS(vsd, file.path(rnaDir, "Metastatic_Tumors_VST_data.rds"))
