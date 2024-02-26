## Description: Test association of ecDNA status and clinical features
## Make figures 1d, 1e, 1f (Figure 2a-c in preprint)
## Also have code for Table 1 p-values

library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(lme4)


# Set params
figureDir <- file.path("figures")
dotSize <- 2.5
conditionColors <- c("ecDNA-" = "#BEAED4", "ecDNA+" = "#7FC97F")


# Load data
sampleInfo <- read_csv("data/Sample_Metadata.csv")
ampliconInfo <- read_csv("data/AC_Results_v2.csv")
clinicalInfo <- read_csv("data/Metastatic_Tumor_Clinical_Data.csv")
tmbInfo <- read_excel("data/Thind_WGS_cSCC_TMB.xlsx")

# Load COSMIC Signature data
sigdf <- read_excel("data/Thind_WGS_Cosmic_Signatures.xls")
sigdf <- as.data.frame(sigdf)

rownames(sigdf) <- sigdf$Signatures
sigdf$Signatures <- NULL
sigdf <- as.matrix(sigdf)
sigdf <- t(sigdf)
sigdf <- as.data.frame(sigdf)

sigdf <- sigdf %>%
  tibble::rownames_to_column("SampleID") %>%
  select(
    SampleID,
    everything()
  )


# Get samples with ecDNA amplicons
ecdnaSamples <- ampliconInfo %>%
  filter(Classification == "ecDNA") %>%
  pull(`Sample name`)

# Label samples as ecDNA- vs ecDNA+
metPatientInfo <- sampleInfo %>%
  inner_join(clinicalInfo, by = c("PatientID" = "Patient")) %>%
  filter(sample_type == "Metastatic Tumor") %>%
  mutate(ecDNA_Status = case_when(
    SampleName %in% ecdnaSamples ~ "ecDNA+",
    TRUE ~ "ecDNA-"
  )) %>%
  inner_join(
    tmbInfo,
    by = c("SampleName" = "Sample ID")
  ) %>%
  inner_join(
    sigdf,
    by = c("SampleName" = "SampleID")
  )

# Create lymph node ratio columns
s <- str_split(metPatientInfo$Lymph_Node_Ratio, "/", simplify = T)
metPatientInfo$Positive_LN <- as.numeric(s[, 1])
metPatientInfo$Total_LN <- as.numeric(s[, 2])
metPatientInfo$LN_Ratio <- metPatientInfo$Positive_LN / metPatientInfo$Total_LN

# Primary and metastatic location have too many categories - skip testing these
# Nodal stage also has many categories but still test

# Test association with age
t.test(Age_Years ~ ecDNA_Status, data=metPatientInfo)
boxplot(Age_Years ~ ecDNA_Status, data=metPatientInfo)

# Test association with Nodal stage
# Use original detailed staging
table(metPatientInfo$Nodal_Stage, metPatientInfo$ecDNA_Status)
fisher.test(table(metPatientInfo$Nodal_Stage, metPatientInfo$ecDNA_Status))

# Collapse nodal staging
metPatientInfo <- metPatientInfo %>%
  mutate(collapsedNodalStage = case_when(
    Nodal_Stage %in% c("N2a", "N2b", "N2c") ~ "N2",
    Nodal_Stage == "N3b" ~ "N3",
    TRUE ~ Nodal_Stage
  ))
table(metPatientInfo$collapsedNodalStage, metPatientInfo$ecDNA_Status)
fisher.test(table(metPatientInfo$collapsedNodalStage, metPatientInfo$ecDNA_Status))


# Test association with extracapsular spread
table(metPatientInfo$Extracapsular_Spread, metPatientInfo$ecDNA_Status)
fisher.test(table(metPatientInfo$Extracapsular_Spread, metPatientInfo$ecDNA_Status))
#
# Remove missing values
with(
  metPatientInfo %>% filter(Extracapsular_Spread != "Not stated"),
  # table(Extracapsular_Spread, ecDNA_Status)
  fisher.test(table(Extracapsular_Spread, ecDNA_Status))
)

# Test association with Grade
table(metPatientInfo$Grade, metPatientInfo$ecDNA_Status)
fisher.test(table(metPatientInfo$Grade, metPatientInfo$ecDNA_Status))

# Remove missing values
with(
  metPatientInfo %>% filter(Grade != "Not stated"),
  # table(Grade, ecDNA_Status)
  fisher.test(table(Grade, ecDNA_Status))
)

# Test association with immunosuppression
table(metPatientInfo$Immunosuppression, metPatientInfo$ecDNA_Status)
fisher.test(table(metPatientInfo$Immunosuppression, metPatientInfo$ecDNA_Status))

# Collapse IS categories
metPatientInfo <- metPatientInfo %>%
  mutate(Immunosuppression_Short = case_when(
    Immunosuppression %in% c("Azathioprine", "cyclosporine A ,tacrolimus") ~ "IS",
    Immunosuppression == "no" ~ "IC"
  ))
table(metPatientInfo$Immunosuppression_Short, metPatientInfo$ecDNA_Status)
fisher.test(table(metPatientInfo$Immunosuppression_Short, metPatientInfo$ecDNA_Status))

# Test association with Lymph Node Ratio
t.test(LN_Ratio ~ ecDNA_Status, data = metPatientInfo)
boxplot(LN_Ratio ~ ecDNA_Status, data = metPatientInfo)

# Better to use aggregated LR model - can take advantage of 
# count info with positive/negative number of lymph nodes
# See https://stats.stackexchange.com/questions/634885/beta-regression-with-success-and-failure-raw-data/634908#634908

respmat <- cbind(metPatientInfo$Positive_LN, metPatientInfo$Total_LN - metPatientInfo$Positive_LN)
f <- glm(respmat ~ ecDNA_Status, data = metPatientInfo, family = binomial)
summary(f)

# Use mixed effects model to account for patient variation
f2 <- glmer(respmat ~ ecDNA_Status + (1 | PatientID), data = metPatientInfo, family = 'binomial')
summary(f2)


## Figure 1d (2a in preprint)
set.seed(52)
lnrPlot <- metPatientInfo %>%
  ggplot(aes(ecDNA_Status, LN_Ratio)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = ecDNA_Status), pch = 21, size = dotSize, position = 'jitter') +
  theme_classic() +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Lymph Node Ratio") +
  theme(
    text = element_text(size = 14)
  ) +
  scale_fill_manual(values = conditionColors)
lnrPlot

## Figure 1e (Figure 2b in preprint)
# Test association with TMB
tmbPlot <- metPatientInfo %>%
  ggplot(aes(ecDNA_Status, `SNVs per Megabase`)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = ecDNA_Status), pch = 21, size = dotSize, position = 'jitter')  +
  theme_classic() +
  labs(x = "", y = "Tumor Mutational Burden (SNVs/Mb)") +
  theme(text = element_text(size = 14)) +
  guides(fill = "none", color = "none") +
  scale_fill_manual(values = conditionColors)
tmbPlot
t.test(`SNVs per Megabase` ~ ecDNA_Status, data = metPatientInfo)

# Test association with UV mutational signatures
# SBS7a-d are UV light exposure. SBS32 is azathioprine
# Note the only sample with SBS32 influence is the patient on azathioprine
boxplot(SBS32 ~ ecDNA_Status, data = metPatientInfo)
t.test(SBS32 ~ ecDNA_Status, data = metPatientInfo)

boxplot(SBS7a + SBS7b + SBS7c + SBS7d ~ ecDNA_Status, data = metPatientInfo)
t.test(SBS7a + SBS7b + SBS7c + SBS7d ~ ecDNA_Status, data = metPatientInfo)

## Figure 1f (Figure 2c in preprint)
metPatientInfo <- metPatientInfo %>%
  mutate(SBS7 = SBS7a + SBS7b + SBS7c + SBS7d)
set.seed(39)
sbs7ScorePlot <- metPatientInfo %>%
  mutate(SBS7 = SBS7a + SBS7b + SBS7c + SBS7d) %>%
  ggplot(aes(ecDNA_Status, SBS7)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = ecDNA_Status), pch = 21, size = dotSize, position = 'jitter') +
  theme_classic() +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "SBS7 Signature") +
  theme(
    text = element_text(size = 14)
  ) +
  scale_fill_manual(values = conditionColors)
sbs7ScorePlot
t.test(SBS7 ~ ecDNA_Status, data = metPatientInfo)

comboPlot <- lnrPlot | tmbPlot | sbs7ScorePlot + plot_layout(guides = "collect")

## Save figures
ggsave(
  filename = file.path(figureDir, "Figure2_LNR_Plot.svg"),
  plot = lnrPlot,
  width = 6,
  height = 8
)
ggsave(
  filename = file.path(figureDir, "Figure2_TMB_Plot.svg"),
  plot = tmbPlot, 
  width = 6, 
  height = 8
)
ggsave(
  filename = file.path(figureDir, "Figure2_SBS7_Plot.svg"),
  plot = sbs7ScorePlot, 
  width = 6, 
  height = 8
)
ggsave(
  filename = file.path(figureDir, "Figure2_Combined.svg"),
  plot = comboPlot,
  width = 16,
  height = 6
)

