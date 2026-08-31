# AWI-Gen phenotype and QC preparation

This script prepares the AWI-Gen phenotype dataset for downstream HIV-1 GWAS analyses.

The workflow:

1. Reads the list of related individuals identified for removal.
2. Reads the AWI-Gen ancestry/admixture results.
3. Reads the AWI-Gen phenotype data.
4. Merges ancestry and phenotype information.
5. Assigns study sites to geographic regions.
6. Recodes HIV status for the Nanoro and Navrongo sites.
7. Applies the ancestry, relatedness and HIV-status QC criteria.
8. Creates the final phenotype files for all participants and separately by sex.

```r
# ============================================================
# AWI-Gen HIV-1 GWAS
# Phenotype and QC preparation
# ============================================================

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

# List of related individuals identified for removal.
# Individuals present in this file are excluded from the
# downstream phenotype dataset.
Dataking <- read.table(
    "/spaces/jeantristan/recall_awigen_ark_jcs/awigen/king/relat_2_arrayunrelated_toberemoved.txt"
)

# AWI-Gen admixture/ancestry results.
DataAdmix <- read.csv(
    "/spaces/jeantristan/recall_awigen_ark_jcs/awigen/admixture/awigen_qcgenet_v3/plot/combineres.tsv",
    sep = "\t"
)

# Alternative admixture result file used previously.
# DataAdmix <- read.csv(
#     "/spaces/jeantristan/recall_awigen_ark_jcs/awigen/admixture/awigen_qcgenet_v3/plot/cleanref_resumealltest.tsv",
#     sep = "\t"
# )

# Raw AWI-Gen phenotype data.
DataPheno <- read.csv(
    "/dataE/AWIGenGWAS/pheno/raw/pheno_2020_07_30.csv"
)


# ------------------------------------------------------------
# Merge ancestry and phenotype information
# ------------------------------------------------------------

# Merge admixture and phenotype data using the study ID.
# Only individuals present in both datasets are retained.
Data2 <- merge(
    DataAdmix,
    DataPheno,
    by.x = "FID",
    by.y = "study_id",
    all = FALSE
)


# ------------------------------------------------------------
# Assign study sites and geographic regions
# ------------------------------------------------------------

# Convert the numerical site code into the corresponding
# AWI-Gen study site.
Data2$site_qc <- c(
    "Agincourt",
    "Dikgale",
    "Nairobi",
    "Nanoro",
    "Navrongo",
    "Soweto"
)[Data2$site]

# Initialise the geographic region.
Data2$region_qc <- "None"

# South Africa sites
Data2$region_qc[
    Data2$site_qc %in% c(
        "Agincourt",
        "Dikgale",
        "Soweto"
    )
] <- "South"

# East Africa
Data2$region_qc[
    Data2$site_qc %in% c("Nairobi")
] <- "East"

# West Africa
Data2$region_qc[
    Data2$site_qc %in% c(
        "Navrongo",
        "Nanoro"
    )
] <- "West"


# ------------------------------------------------------------
# HIV phenotype recoding
# ------------------------------------------------------------

# Recode HIV status for participants from Nanoro and Navrongo.
# Participants with hiv_final_status_c == 2 are recoded as 0.
Data2$hiv_final_status_c[
    Data2$site_qc %in% c("Nanoro", "Navrongo") &
    Data2$hiv_final_status_c == 2
] <- 0


# ------------------------------------------------------------
# Ancestry QC
# ------------------------------------------------------------

# Flag Nairobi participants with an African ancestry proportion
# greater than 0.70.
#
# These individuals are marked in NoAdm and are retained by the
# final QC filter below.
Data2$NoAdm[
    Data2$site_qc == "Nairobi" &
    Data2$Adm_AFR > 0.7
] <- TRUE


# ------------------------------------------------------------
# Sex
# ------------------------------------------------------------

# Sex coding:
# 0 = Female
# 1 = Male
Data2$sex_qc <- NA_character_

Data2$sex_qc[Data2$sex == 0] <- "Female"
Data2$sex_qc[Data2$sex == 1] <- "Male"


# ------------------------------------------------------------
# Final phenotype/QC filter
# ------------------------------------------------------------

# Retain individuals who:
#   1. Have available EigPcsFilter information;
#   2. Are not listed among related individuals to be removed;
#   3. Do not have HIV status coded as 2;
#   4. Have a TRUE value for NoAdm.
#
# This defines the final dataset used for downstream GWAS.
balise1 <- (
    !is.na(Data2$EigPcsFilter) &
    !(Data2$FID %in% Dataking$V1) &
    Data2$hiv_final_status_c != 2 &
    (!is.na(Data2$NoAdm) & Data2$NoAdm)
)

Data3 <- Data2[balise1, ]


# ------------------------------------------------------------
# Final GWAS variables
# ------------------------------------------------------------

# Rename the final HIV status variable used for GWAS.
Data3$hiv <- Data3$hiv_final_status_c

# Generate an age-squared term.
# This can be used to model a potential non-linear relationship
# between age and HIV status in downstream analyses.
Data3$age2 <- Data3$age ** 2


# ------------------------------------------------------------
# Generate phenotype files
# ------------------------------------------------------------

# Generate one file for each geographic region and additional
# sex-specific files within each region.
#
# Regional files:
#   South.all
#   East.all
#   West.all
#
# Sex-specific files:
#   South_Female.gender
#   South_Male.gender
#   East_Female.gender
#   East_Male.gender
#   West_Female.gender
#   West_Male.gender

for (region in unique(Data3$region_qc)) {

    # All participants within the region
    write.table(
        Data3[
            Data3$region_qc == region,
            c(
                "FID",
                "IID",
                "hiv",
                "age",
                "age2",
                "sex_qc"
            )
        ],
        file = paste(region, ".all", sep = ""),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )

    # Female and male participants separately
    for (gender in unique(Data3$sex_qc)) {

        write.table(
            Data3[
                Data3$sex_qc == gender &
                Data3$region_qc == region,
                c(
                    "FID",
                    "IID",
                    "hiv",
                    "age",
                    "age2",
                    "sex_qc"
                )
            ],
            file = paste(
                region,
                "_",
                gender,
                ".gender",
                sep = ""
            ),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
    }
}
```

## Output

The script generates phenotype files by geographic region and sex, including:

```text
South.all
East.all
West.all

South_Female.gender
South_Male.gender

East_Female.gender
East_Male.gender

West_Female.gender
West_Male.gender
```

These files provide the phenotype and covariate information used in the downstream AWI-Gen HIV-1 GWAS analyses.

