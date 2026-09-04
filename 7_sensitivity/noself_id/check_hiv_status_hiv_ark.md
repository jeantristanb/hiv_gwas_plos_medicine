---
title: "Prepare ARK HIV Phenotype Files"
author: "Jean-Tristan Brandenburg"
date: "2026-08-05"
output: html_document
---

# Overview

This document prepares HIV phenotype files for downstream genetic analyses.

The workflow:

1. Loads ARK phenotype data.
2. Loads ancestry/admixture quality-control results.
3. Selects individuals passing ancestry QC.
4. Restricts to unrelated individuals.
5. Removes participants with missing or unknown HIV status.
6. Creates standardized phenotype files.
7. Exports overall, female-only, and male-only datasets.


``` r
library(dplyr)
```

# Input data


``` r
# ARK phenotype data
Data1 <- read.csv(
  "/dataA/ark/pheno/Phase1ARKStudyAginco_DATA_2024-01-30_1333.csv"
)

# Admixture/QC results
DataAdmixture <- read.csv(
  "/spaces/jeantristan/recall_awigen_ark_jcs/ark/admixture/combineres.tsv",
  sep = "\t"
)

# List of unrelated individuals
unrelated <- read.table(
  "/spaces/jeantristan/recall_awigen_ark_jcs/all/king/all_unrelated_2unrelated.txt"
)
```

# Select individuals passing ancestry QC

Individuals are retained if they:

- passed the PCA filter (`PassPcsFilt == TRUE`)
- were classified as non-admixed (`NoAdm == TRUE`)


``` r
keep_admixture <-
  !is.na(DataAdmixture$PassPcsFilt) &
  DataAdmixture$PassPcsFilt &
  DataAdmixture$NoAdm

sample_adm <- DataAdmixture$FID[keep_admixture]
```

# Select eligible participants

Participants are retained if they:

- have a known HIV status (`arkp1_hiv_status != 2`)
- passed ancestry QC
- are unrelated according to KING
- [arkp1_hiv_status] 	Based on (i) prior HIV testing history OR (ii) HIV PCR testing for ARK this participants HIV status is:	negative & positive


``` r
keep <- !is.na(Data1$arkp1_hiv_status) &
  Data1$arkp1_hiv_status != 2 &
  Data1$arkp1bcode %in% sample_adm &
  Data1$arkp1bcode %in% unrelated$V1
```

# Create phenotype variables

Sex coding is converted to PLINK convention:

| Original | New |
|----------:|----:|
| Male (1) | 1 |
| Female (2) | 0 |


``` r
Data1$sex <- Data1$arkp1_sex
Data1$sex[Data1$sex == 2] <- 0

Data1$age <- Data1$arkp1_age
```

# Create phenotype table


``` r
columns <- c(
  "arkp1bcode",
  "arkp1bcode",
  "arkp1_hiv_status",
  "sex",
  "arkp1_age"
)

new_names <- c(
  "FID",
  "IID",
  "hiv",
  "sex",
  "age"
)

DataF <- Data1[keep, columns]
names(DataF) <- new_names
```

# Summary


``` r
cat("Total participants:", nrow(DataF), "\n")
```

```
## Total participants: 654
```

``` r
cat("Females:", sum(DataF$sex == 0), "\n")
```

```
## Females: 402
```

``` r
cat("Males:", sum(DataF$sex == 1), "\n")
```

```
## Males: 252
```

``` r
table(DataF$hiv)
```

```
## 
##   0   1 
## 515 139
```

# Export phenotype files


``` r
#write.table(DataF, file = "ark_hiv.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE,sep = "\t")

#write.table(
#  DataF[DataF$sex == 0, ],
#  file = "ark_hiv_females.tsv",
#  row.names = FALSE,
#  col.names = TRUE,
#  quote = FALSE,
#  sep = "\t"
#)

#write.table(
#  DataF[DataF$sex == 1, ],
#  file = "ark_hiv_males.tsv",
#  row.names = FALSE,
#  col.names = TRUE,
#  quote = FALSE,
#  sep = "\t"
#)
```

## Compare reported HIV status with HIV test results

### Create HIV status variables

The reported HIV status is taken from `arkp1_hiv_dx`, where values greater than 1 are considered missing.


``` r
# Reported HIV status
Data1$hiv_positive <- Data1$arkp1_hiv_dx
Data1$hiv_positive[Data1$hiv_positive > 1] <- NA

# HIV test result
Data1$result_hiv <- NA

# At least one HIV test available
has_result <- !is.na(Data1$arkp1_hiv1_res) |
              !is.na(Data1$arkp1_hiv2_res)

Data1$result_hiv[has_result] <- 0

# Positive if both HIV tests are positive
positive_result <-
  !is.na(Data1$arkp1_hiv1_res) &
  Data1$arkp1_hiv1_res == 1 &
  !is.na(Data1$arkp1_hiv2_res) &
  Data1$arkp1_hiv2_res == 1

Data1$result_hiv[positive_result] <- 1
```

### Compare reported HIV status with laboratory results


``` r
library(dplyr)
library(tidyr)
library(knitr)

Data1 <- Data1[Data1$arkp1bcode %in% DataF$FID, ]

comparison <- Data1 %>%
  mutate(
    result_hiv = factor(
      result_hiv,
      levels = c(0, 1, NA),
      labels = c("Negative", "Positive", "Missing"),
      exclude = NULL
    )
  ) %>%
  group_by(hiv_positive, result_hiv) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = result_hiv,
    values_from = N,
    values_fill = 0
  ) %>%
  rename(
    Initial_status = hiv_positive
  ) %>%
  mutate(
    Total = Negative + Positive + Missing,
    `Percent Missing` = round(100 * Missing / Total, 1)
  )

# Add total row
comparison <- comparison %>%
  mutate(Initial_status = as.character(Initial_status))

comparison <- bind_rows(
  comparison,
  comparison %>%
    summarise(
      Initial_status = "Total",
      Negative = sum(Negative),
      Positive = sum(Positive),
      Missing = sum(Missing),
      Total = sum(Total),
      `Percent Missing` = round(100 * sum(Missing) / sum(Total), 1)
    )
)

kable(comparison)
```



|Initial_status | Negative| Positive| Missing| Total| Percent Missing|
|:--------------|--------:|--------:|-------:|-----:|---------------:|
|0              |      135|        6|     360|   501|            71.9|
|1              |        0|        0|     131|   131|           100.0|
|NA             |       20|        2|       0|    22|             0.0|
|Total          |      155|        8|     491|   654|            75.1|

``` r
write.csv(comparison,file='comp_ark.csv',row.names=F)
```

Retain only individuals with an available HIV test result (result_hiv). Exclude individuals with an initial HIV-positive status (hiv_positive = 1) but a negative study HIV test result (result_hiv = 0).


``` r
# Retain participants with an available HIV test result and remove
# discordant individuals with an initial HIV-positive status but a
# negative HIV test result.

Data_F_sel <- Data1 %>%
  filter(
    !is.na(result_hiv),
    !(hiv_positive == 1 & result_hiv == 0)
  )

columns <- c(
  "arkp1bcode",
  "arkp1bcode",
  "arkp1_hiv_status",
  "sex",
  "arkp1_age"
)

new_names <- c(
  "FID",
  "IID",
  "hiv",
  "sex",
  "age"
)

Data_F_sel<- Data_F_sel[, columns]
names(Data_F_sel) <- new_names
dir.create('ark_qc_justtest/')
```

```
## Warning in dir.create("ark_qc_justtest/"): 'ark_qc_justtest' already exists
```

``` r
write.table(
  Data_F_sel[Data_F_sel$sex == 0, ],
  file = "ark_qc_justtest/ark_hiv_females.tsv",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)

write.table(
  Data_F_sel[Data_F_sel$sex == 1, ],
  file = "ark_qc_justtest/ark_hiv_males.tsv",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)
```
