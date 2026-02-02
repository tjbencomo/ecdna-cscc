## Code for analysis of ecDNA associations with clinical features

library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)

########################################################
####### Set params
########################################################
figureDir <- file.path("figures")
dotSize <- 2.5
conditionColors <- c("ecDNA-" = "#BEAED4", "ecDNA+" = "#7FC97F")

########################################################
####### Load data
########################################################
sampleInfo <- read_csv("data/Sample_Metadata.csv")
ampliconInfo <- read_csv("data/AC_Results_v2.csv")
clinicalInfo <- read_csv("data/Metastatic_Tumor_Clinical_Data.csv")
tmbInfo <- read_excel("data/Thind_WGS_cSCC_TMB.xlsx")

pni_lvi_info <- read_excel("data/UOW_BA_Tumor_Data.xlsx")
pni_lvi_info <- pni_lvi_info %>%
  select(CSCC, `M/N PNI`, `M/N LVI`) %>%
  rename(
    Patient = CSCC,
    PNI_Status = `M/N PNI`,
    LVI_Status = `M/N LVI`
  )

# Remove CSCC_0024 because wasn't run with AA due to errors
clinicalInfo <- clinicalInfo %>%
  filter(Patient != "CSCC_0024") %>%
  left_join(pni_lvi_info)


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

########################################################
####### Associations for Table 2
########################################################

# Test association with Nodal stage
metPatientInfo <- metPatientInfo %>%
  mutate(collapsedNodalStage = case_when(
    Nodal_Stage %in% c("N2a", "N2b", "N2c") ~ "N2",
    Nodal_Stage %in% c("3b", "N3b") ~ "N3",
    TRUE ~ Nodal_Stage
  ))
table(metPatientInfo$collapsedNodalStage, metPatientInfo$ecDNA_Status)
fisher.test(table(metPatientInfo$collapsedNodalStage, metPatientInfo$ecDNA_Status))


# Test association with extracapsular spread
with(
  metPatientInfo %>% filter(Extracapsular_Spread != "Not stated"),
  {
    print(table(Extracapsular_Spread, ecDNA_Status))
    print(fisher.test(table(Extracapsular_Spread, ecDNA_Status)))
  }
)

# Test association with Grade
with(
  metPatientInfo %>% filter(Grade != "Not stated"),
  {
    print(table(Grade, ecDNA_Status))
    print(fisher.test(table(Grade, ecDNA_Status)))
  }
)

# Test association with immunosuppression
metPatientInfo <- metPatientInfo %>%
  mutate(Immunosuppression_Short = case_when(
    Immunosuppression %in% c("Azathioprine", "cyclosporine A ,tacrolimus") ~ "IS",
    Immunosuppression == "no" ~ "IC"
  ))
table(metPatientInfo$Immunosuppression_Short, metPatientInfo$ecDNA_Status)
fisher.test(table(metPatientInfo$Immunosuppression_Short, metPatientInfo$ecDNA_Status))

# Test association with PNI
table(metPatientInfo$PNI_Status, metPatientInfo$ecDNA_Status)
fisher.test(table(metPatientInfo$PNI_Status, metPatientInfo$ecDNA_Status))

# Test association with LVI
table(metPatientInfo$LVI_Status, metPatientInfo$ecDNA_Status)
fisher.test(table(metPatientInfo$LVI_Status, metPatientInfo$ecDNA_Status))


########################################################
####### Figure 1C and 1G
########################################################

# Test association with TMB
figure1C_1 <- metPatientInfo %>%
  ggplot(aes(ecDNA_Status, `SNVs per Megabase`)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = ecDNA_Status), pch = 21, size = dotSize, position = 'jitter')  +
  theme_classic() +
  labs(x = "", y = "Tumor Mutational Burden (SNVs/Mb)") +
  theme(text = element_text(size = 14)) +
  guides(fill = "none", color = "none") +
  scale_fill_manual(values = conditionColors)
figure1C_1
t.test(`SNVs per Megabase` ~ ecDNA_Status, data = metPatientInfo)

# Test association with UV mutational signatures
# SBS7a-d are UV light exposure. SBS32 is azathioprine
# Note the only sample with SBS32 influence is the patient on azathioprine
metPatientInfo <- metPatientInfo %>%
  mutate(SBS7 = SBS7a + SBS7b + SBS7c + SBS7d)

set.seed(39)
figure1C_2 <- metPatientInfo %>%
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
figure1C_2
t.test(SBS7 ~ ecDNA_Status, data = metPatientInfo)


# Test association with Lymph Node Ratio
# Better to use aggregated LR model - can take advantage of 
# count info with positive/negative number of lymph nodes
# See https://stats.stackexchange.com/questions/634885/beta-regression-with-success-and-failure-raw-data/634908#634908
respmat <- cbind(metPatientInfo$Positive_LN, metPatientInfo$Total_LN - metPatientInfo$Positive_LN)
f <- glm(respmat ~ ecDNA_Status, data = metPatientInfo, family = quasibinomial())
summary(f)

set.seed(52)
figure1G <- metPatientInfo %>%
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
figure1G

# Test association with total number of lymph nodes tested
f_total <- lm(Total_LN ~ ecDNA_Status, data = metPatientInfo)
summary(f_total)