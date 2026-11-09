---
title: "AWI-Gen Phenotype Processing for HIV Analysis"
author: "Jean-Tristan Brandenburg"
date: "2026-08-05"
output:
  html_document:
    toc: true
    toc_depth: 3
    number_sections: true
  pdf_document: default
---

# Overview

This document prepares the AWI-Gen phenotype dataset for analyses involving HIV status.

The workflow includes:

- Loading phenotype and admixture data
- Merging datasets
- Annotating study sites and geographical regions
- Harmonizing HIV status
- Applying genetic quality-control filters
- Creating derived variables
- Recovering original HIV variables
- Excluding West African participants
- Cleaning original HIV variables


``` r
library(dplyr)
library(tidyr)
```

# Load phenotype datasets


``` r
# Processed phenotype dataset
DataPheno <- read.csv(
  "/dataE/AWIGenGWAS/pheno/raw/pheno_2020_07_30.csv"
)

Dataking<-read.table('/spaces/jeantristan/recall_awigen_ark_jcs/awigen//king/relat_2_arrayunrelated_toberemoved.txt')
DataAdmix<-read.csv('/spaces/jeantristan/recall_awigen_ark_jcs/awigen/admixture/awigen_qcgenet_v3/plot/combineres.tsv',sep='\t')


# Original phenotype database
newdata_pheno <- read.csv(
  "/spaces/jeantristan/ckd_project/DataI/all_awigen1.csv"
)
```

# Merge admixture and phenotype data

The phenotype information is merged with the genetic dataset using the participant identifier.


``` r
Data2 <- merge(
  DataAdmix,
  DataPheno,
  by.x = "FID",
  by.y = "study_id",
  all = FALSE
)
```

# Annotate study site and region

Study sites are converted from numeric codes into readable names and grouped into geographical regions.


``` r
Data2$site_qc <- c(
  "Agincourt",
  "Dikgale",
  "Nairobi",
  "Nanoro",
  "Navrongo",
  "Soweto"
)[Data2$site]

Data2$region_qc <- "None"

Data2$region_qc[
  Data2$site_qc %in% c("Agincourt", "Dikgale", "Soweto")
] <- "South"

Data2$region_qc[
  Data2$site_qc == "Nairobi"
] <- "East"

Data2$region_qc[
  Data2$site_qc %in% c("Nanoro", "Navrongo")
] <- "West"

table(Data2$site_qc)
```

```
## 
## Agincourt   Dikgale   Nairobi    Nanoro  Navrongo    Soweto 
##      2304      1141      1952      1992      1858      1761
```

``` r
table(Data2$region_qc)
```

```
## 
##  East South  West 
##  1952  5206  3850
```

# HIV status harmonization

The variable **hiv_final_status_c** summarizes HIV status.

| Value | Interpretation |
|------:|----------------|
|0|HIV negative|
|1|HIV positive|
|2|Unknown|

The original derivation was

```text
CASE
WHEN (result_hiv IS NULL OR result_hiv = 2)
 AND (hiv_positive IS NULL OR hiv_positive >=2)
THEN 2

WHEN result_hiv = 1
 OR hiv_positive = 1
THEN 1

ELSE 0
END
```

For Nanoro and Navrongo, unknown HIV status was treated as HIV negative.


``` r
Data2$hiv_final_status_c[
  Data2$site_qc %in% c("Nanoro", "Navrongo") &
  Data2$hiv_final_status_c == 2
] <- 0

table(Data2$hiv_final_status_c)
```

```
## 
##    0    1    2 
## 8943 1381  684
```

# Ancestry and sex variables

Only participants with predominantly African ancestry are retained for the Nairobi cohort.


``` r
Data2$NoAdm[
  Data2$site_qc == "Nairobi" &
  Data2$Adm_AFR > 0.7
] <- TRUE

Data2$sex_qc <- NA_character_

Data2$sex_qc[Data2$sex == 0] <- "Female"
Data2$sex_qc[Data2$sex == 1] <- "Male"

table(Data2$sex_qc)
```

```
## 
## Female   Male 
##   6026   4982
```

# Apply quality-control filters

Participants are retained if they satisfy all of the following:

- Pass principal component filtering
- Are not duplicated/related (`Dataking`)
- Have known HIV status
- Pass ancestry filtering


``` r
keep <- !is.na(Data2$EigPcsFilter) &
  !(Data2$FID %in% Dataking$V1) &
  Data2$hiv_final_status_c != 2 &
  !is.na(Data2$NoAdm) &
  Data2$NoAdm

Data3 <- Data2[keep, ]

dim(Data3)
```

```
## [1] 8906  187
```

# Create derived variables


``` r
Data3$hiv <- Data3$hiv_final_status_c

Data3$age2 <- Data3$age^2

summary(Data3$age)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   29.00   45.00   51.00   51.55   56.00   80.00
```

# Recover original HIV variables

The original HIV variables are reintroduced from the raw phenotype database.


``` r
newdata_pheno2 <- newdata_pheno %>%
  select(study_id, result_hiv, hiv_positive)

newdata_pheno2$in_db <- TRUE

Data3F <- merge(
  Data3,
  newdata_pheno2,
  by.x = "FID",
  by.y = "study_id",
  all.x = TRUE
)
```

# Exclude West African participants

For analyses restricted to East and Southern African populations.


``` r
Data3F<- subset(
  Data3F,
  region_qc != "West"
)
```

## Clean original HIV variables

The original phenotype database contains two HIV-related variables:

- **`hiv_positive`**: Initial HIV status recorded at enrolment, based on self-report and/or existing clinical records.
- **`result_hiv`**: HIV test result obtained during the AWI-Gen study.

Both variables use the following coding:

- **0** = HIV negative
- **1** = HIV positive
- **≥2** = Unknown or unavailable HIV status

For all analyses, values **≥2** were recoded to `NA` to represent missing data.


``` r
Data3F$hiv_positive[
  Data3F$hiv_positive >= 2
] <- NA

Data3F$result_hiv[
  Data3F$result_hiv >= 2
] <- NA
```


``` r
hiv_summary <- Data3F %>%
  mutate(
    hiv_positive_available = !is.na(hiv_positive),
    result_hiv_available = !is.na(result_hiv),
    both_available = hiv_positive_available & result_hiv_available
  ) %>%
  group_by(site_qc) %>%
  summarise(
    N = n(),
    HIV_positive_non_missing = sum(hiv_positive_available),
    HIV_positive_percent = round(100 * mean(hiv_positive_available), 1),
    Result_HIV_non_missing = sum(result_hiv_available),
    Result_HIV_percent = round(100 * mean(result_hiv_available), 1),
    Both_non_missing = sum(both_available),
    Both_percent = round(100 * mean(both_available), 1),
    .groups = "drop"
  )


kable(hiv_summary)
```



|site_qc   |    N| HIV_positive_non_missing| HIV_positive_percent| Result_HIV_non_missing| Result_HIV_percent| Both_non_missing| Both_percent|
|:---------|----:|------------------------:|--------------------:|----------------------:|------------------:|----------------:|------------:|
|Agincourt | 1768|                     1765|                 99.8|                   1676|               94.8|             1673|         94.6|
|Dikgale   |  935|                      718|                 76.8|                    784|               83.9|              567|         60.6|
|Nairobi   | 1613|                     1531|                 94.9|                    579|               35.9|              497|         30.8|
|Soweto    | 1224|                      558|                 45.6|                   1048|               85.6|              382|         31.2|


# Comparison between initial status and result hiv


``` r
library(dplyr)
library(tidyr)

comparison <- Data3F %>%
  mutate(
    Initial_status = factor(
      hiv_positive,
      levels = c(0, 1, NA),
      labels = c("Initial negative", "Initial positive", "Initial missing"),
      exclude = NULL
    ),
    result_hiv = factor(
      result_hiv,
      levels = c(0, 1, NA),
      labels = c("Result negative", "Result positive", "Missing"),
      exclude = NULL
    )
  ) %>%
  group_by(site_qc, Initial_status, result_hiv) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = result_hiv,
    values_from = N,
    values_fill = 0
  ) %>%
  mutate(
    Total = `Result negative` + `Result positive` + Missing,
    `% Missing` = round(100 * Missing / Total, 1)
  )

# Add totals for each site
comparison <- comparison %>%
  group_by(site_qc) %>%
  group_modify(~{
    total_row <- .x %>%
      summarise(
        Initial_status = "Total",
        `Result negative` = sum(`Result negative`),
        `Result positive` = sum(`Result positive`),
        Missing = sum(Missing),
        Total = sum(Total),
        `% Missing` = round(100 * Missing / Total, 1)
      )

    bind_rows(.x, total_row)
  }) %>%
  ungroup()

# Overall summary (excluding site totals)
overall <- comparison %>%
  filter(Initial_status != "Total") %>%
  group_by(Initial_status) %>%
  summarise(
    site_qc = "All sites",
    `Result negative` = sum(`Result negative`),
    `Result positive` = sum(`Result positive`),
    Missing = sum(Missing),
    Total = sum(Total),
    `% Missing` = round(100 * Missing / Total, 1),
    .groups = "drop"
  )

# Overall total
overall_total <- overall %>%
  summarise(
    site_qc = "All sites",
    Initial_status = "Total",
    `Result negative` = sum(`Result negative`),
    `Result positive` = sum(`Result positive`),
    Missing = sum(Missing),
    Total = sum(Total),
    `% Missing` = round(100 * Missing / Total, 1)
  )

# Combine everything
comparison <- bind_rows(
  comparison,
  overall,
  overall_total
)

knitr::kable(comparison)
```



|site_qc   |Initial_status   | Result negative| Result positive| Missing| Total| % Missing|
|:---------|:----------------|---------------:|---------------:|-------:|-----:|---------:|
|Agincourt |Initial negative |            1209|             210|      80|  1499|       5.3|
|Agincourt |Initial positive |              13|             241|      12|   266|       4.5|
|Agincourt |Initial missing  |               0|               3|       0|     3|       0.0|
|Agincourt |Total            |            1222|             454|      92|  1768|       5.2|
|Dikgale   |Initial negative |             524|              26|      34|   584|       5.8|
|Dikgale   |Initial positive |               3|              14|     117|   134|      87.3|
|Dikgale   |Initial missing  |             193|              24|       0|   217|       0.0|
|Dikgale   |Total            |             720|              64|     151|   935|      16.1|
|Nairobi   |Initial negative |             463|              29|     849|  1341|      63.3|
|Nairobi   |Initial positive |               1|               4|     185|   190|      97.4|
|Nairobi   |Initial missing  |              77|               5|       0|    82|       0.0|
|Nairobi   |Total            |             541|              38|    1034|  1613|      64.1|
|Soweto    |Initial negative |             357|              15|      57|   429|      13.3|
|Soweto    |Initial positive |               8|               2|     119|   129|      92.2|
|Soweto    |Initial missing  |             552|             114|       0|   666|       0.0|
|Soweto    |Total            |             917|             131|     176|  1224|      14.4|
|All sites |Initial missing  |             822|             146|       0|   968|       0.0|
|All sites |Initial negative |            2553|             280|    1020|  3853|      26.5|
|All sites |Initial positive |              25|             261|     433|   719|      60.2|
|All sites |Total            |            3400|             687|    1453|  5540|      26.2|

``` r
write.csv(comparison, row.names=F, file='comp_awigen.csv')
```

# selection of sample for specific analysis 
 * reviewer comment about self identificat
 * keep just individual with result but exclude  individual with self identificaiton positive but negative 0? -> on medication `
 * Retain only individuals with an available HIV test result (result_hiv). Exclude individuals who self-identified as HIV positive (hiv_positive = 1) but had a negative HIV test result (result_hiv = 0), as these individuals are likely receiving antiretroviral treatment and have a suppressed viral load.
	

``` r
if(!file.exists('pheno_qcnoreported/'))dir.create('pheno_qcnoreported/')
Data_final <- Data3F[ !is.na(Data3F$result_hiv) & (is.na(Data3F$hiv_positive) | !(Data3F$hiv_positive == 1 & Data3F$result_hiv == 0)), ]

for(region in unique(Data_final$region_qc)){
        write.table(Data_final[Data_final$region_qc==region,c('FID','IID','hiv','age','age2','sex_qc')], file=paste('pheno_qcnoreported/',region,'.all',sep=''),sep='\t',row.names=F, col.names=T, quote=F)
        for(gender in unique(Data_final$sex_qc))write.table(Data_final[Data_final$sex_qc==gender & Data_final$region_qc==region,c('FID','IID','hiv','age','age2','sex_qc')], file=paste('pheno_qcnoreported/',region,'_',gender,'.gender',sep=''),sep='\t',row.names=F, col.names=T, quote=F)
}
```



# Session information


``` r
sessionInfo()
```

```
## R version 4.6.1 (2026-06-24)
## Platform: x86_64-redhat-linux-gnu
## Running under: Rocky Linux 9.8 (Blue Onyx)
## 
## Matrix products: default
## BLAS/LAPACK: FlexiBLAS OPENBLAS-OPENMP;  LAPACK version 3.9.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Africa/Johannesburg
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] tidyr_1.3.2 dplyr_1.2.1 knitr_1.51 
## 
## loaded via a namespace (and not attached):
##  [1] R6_2.6.1         tidyselect_1.2.1 xfun_0.60        magrittr_2.0.5  
##  [5] glue_1.8.1       tibble_3.3.1     pkgconfig_2.0.3  generics_0.1.4  
##  [9] lifecycle_1.0.5  cli_3.6.6        vctrs_0.7.3      withr_3.0.2     
## [13] compiler_4.6.1   purrr_1.2.2      tools_4.6.1      evaluate_1.0.5  
## [17] pillar_1.11.1    otel_0.2.0       rlang_1.2.0
```
