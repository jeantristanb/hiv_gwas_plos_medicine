# ARK HIV-1 phenotype and QC preparation

This script prepares the ARK phenotype dataset for downstream HIV-1 GWAS analyses.

HIV status was defined based on **prior HIV testing history and/or HIV PCR testing** available for ARK participants.

The workflow:

1. Reads the ARK phenotype data.
2. Reads the ancestry/admixture QC results.
3. Selects individuals passing the principal-component/admixture QC.
4. Excludes related individuals using the KING unrelatedness list.
5. Excludes individuals with missing or unresolved HIV status.
6. Recodes sex to the format used for the GWAS.
7. Creates overall, female-specific and male-specific phenotype files.

## Phenotype and QC preparation

```r
# ============================================================
# ARK HIV-1 GWAS
# Phenotype and QC preparation
# ============================================================

# ------------------------------------------------------------
# Input phenotype data
# ------------------------------------------------------------

# ARK Phase 1 phenotype data.
#
# Based on:
#   (i) prior HIV testing history and/or
#   (ii) HIV PCR testing
#
# the ARK HIV status is defined in arkp1_hiv_status.
Data1 <- read.csv(
    "/home/jeantristan/Ark/ark/Phenotype/Phase1ARKStudyAginco_DATA_2024-01-30_1333.csv"
)


# ------------------------------------------------------------
# Ancestry/admixture QC
# ------------------------------------------------------------

# AWI-Gen/ARK admixture and PCA QC results.
DataAdmixture <- read.csv(
    "/spaces/jeantristan/recall_awigen_ark_jcs/ark/admixture/combineres.tsv",
    sep = "\t"
)

# Select individuals who:
#   - have a non-missing PassPcsFilt result;
#   - pass the PCA filter;
#   - pass the ancestry/admixture filter (NoAdm).
balise <- (
    !is.na(DataAdmixture$PassPcsFilt) &
    DataAdmixture$PassPcsFilt == TRUE &
    DataAdmixture$NoAdm
)

# IDs of individuals passing the ancestry/PCA QC.
sample_adm <- DataAdmixture$FID[balise]


# ------------------------------------------------------------
# Relatedness QC
# ------------------------------------------------------------

# List of unrelated individuals identified using KING.
unrelated <- read.table(
    "/spaces/jeantristan/recall_awigen_ark_jcs/all/king/all_unrelated_2unrelated.txt"
)


# ------------------------------------------------------------
# Final participant selection
# ------------------------------------------------------------

# Retain individuals who:
#   1. have a non-missing HIV status;
#   2. have a resolved HIV status (status != 2);
#   3. passed the ancestry/PCA QC;
#   4. are included in the unrelated-individual list.
baliseall <- (
    !is.na(Data1$arkp1_hiv_status) &
    Data1$arkp1_hiv_status != 2 &
    Data1$arkp1bcode %in% sample_adm &
    Data1$arkp1bcode %in% unrelated$V1
)


# ------------------------------------------------------------
# Sex coding
# ------------------------------------------------------------

# Original ARK coding:
#   1 = Male
#   2 = Female
#
# GWAS coding:
#   0 = Female
#   1 = Male

Data1$sex <- Data1$arkp1_sex
Data1$sex[Data1$sex == 2] <- 0


# ------------------------------------------------------------
# Age
# ------------------------------------------------------------

# Use the ARK age variable as the GWAS age covariate.
Data1$age <- Data1$arkp1_age


# ------------------------------------------------------------
# Select and rename variables
# ------------------------------------------------------------

# Variables retained for downstream GWAS:
#
#   FID  = participant identifier
#   IID  = participant identifier
#   hiv  = binary HIV status
#   sex  = sex
#   age  = age

listhead <- c(
    "arkp1bcode",
    "arkp1bcode",
    "arkp1_hiv_status",
    "sex",
    "arkp1_age"
)

listnewhead <- c(
    "FID",
    "IID",
    "hiv",
    "sex",
    "age"
)

DataF <- Data1[
    baliseall,
    listhead
]

names(DataF) <- listnewhead


# ------------------------------------------------------------
# Write final phenotype files
# ------------------------------------------------------------

# All eligible ARK participants
write.table(
    DataF,
    file = "ark_hiv.tsv",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
)

# Female participants
write.table(
    DataF[DataF$sex == 0, ],
    file = "ark_hiv_females.tsv",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
)

# Male participants
write.table(
    DataF[DataF$sex == 1, ],
    file = "ark_hiv_males.tsv",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
)
```

# build pcs 
```sh
/home/jeantristan/bin/plink --bfile ${base} --indep-pairwise 50 5 0.2 --keep $file --out $headout"_"check --maf 0.01
/home/jeantristan/bin/plink --bfile ${base} --extract $headout"_"check.prune.in --make-bed --out $headout"_"check  --maf 0.01 --keep $file
/home/jeantristan/bin/plink --bfile $headout"_"check --pca --out $headout
~/Travail/git/gwas_ckdawigen/bin/addpheno_pcs.r --data $file --pcs $headout".eigenvec" --out $dirpheno/$head".ind"
```

## Output

The script generates three phenotype files:

```text
ark_hiv.tsv
ark_hiv_females.tsv
ark_hiv_males.tsv
```

The female and male files are subsequently used for the sex-specific ARK HIV-1 GWAS analyses.

## QC criteria

The final dataset includes participants who:

* passed the PCA/admixture QC;
* passed the ancestry/admixture criterion;
* were classified as unrelated according to KING;
* had a non-missing HIV status;
* did not have HIV status coded as `2`.

## Sex coding

```text
0 = Female
1 = Male
```

Individual-level ARK phenotype and genotype data are not included in the public repository.

