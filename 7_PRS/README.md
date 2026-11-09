# HIV Polygenic Risk Score Analysis

This repository contains the code used to calculate and evaluate polygenic risk scores (PRS) for HIV acquisition across multiple study sites.

The workflow includes:

1. PRS calculation
2. Principal component (PC) calculation
3. PRS extraction
4. Statistical analysis
5. Validation

The PRS workflow used by `run_prs.bash` is maintained in the GitHub repository **`jeantristanb/qc_gwas_genet`** and is available locally under:

```text
~/Cancer/qc_gwas_genet
```

> **Data access:** Individual-level genotype and phenotype data are not included in this repository. Users must have appropriate authorization to access and analyse the underlying study data.

---

## Repository structure

```text
.
├── README.md
└── code_source/
    ├── build_pcs.r
    ├── checkval.r
    ├── extract_genetics_data.bash
    ├── extract_stat.r
    └── run_prs.bash
```

---

# Analysis workflow

```text
Genotype data
      |
      v
PRS calculation
      |
      v
PRS extraction
      |
      +------------------+
      |                  |
      v                  v
Principal Components   Phenotype data
      |                  |
      +--------+---------+
               |
               v
        Statistical analysis
          /             \
         v               v
 Logistic regression   ROC / AUC
          \             /
           v           v
              Results
                 |
                 v
             Validation
```

---

# Code source

## `run_prs.bash`

Main script used to calculate the HIV PRS.

The workflow uses the HIV GWAS summary statistics together with genotype data and the PRS-CS EUR LD reference.

### Main inputs

- PLINK genotype dataset
- HIV phenotype data
- HIV GWAS summary statistics
- PRS-CS EUR LD reference

### Main analysis options

- MAF threshold: `0.01`
- PRS-CS enabled
- Five PCs added by the PRS workflow
- Covariates: `age`, `sex`, `Study.Site`
- Nextflow
- SLURM/Singularity execution

The original workflow uses:

```bash
module load nextflow

refprsc_eur="/path/to/PRScs/ldblk_1kg_eur/"
bfile="/path/to/genotype/all_clean"
dataI="../data/pheno/allsite.tsv"
sumstateur="../external_data/HIV_acquisition_meta_sample_size1_filtered_annotated.txt.hg38.bcftools"
pheno="hiv"
Neur=594646

data=$dataI

~/Cancer/qc_gwas_genet/prs.nf \
    --ref_prscs $refprsc_eur \
    --bfile $bfile \
    --data $data \
    --pheno $pheno \
    --sumstat $sumstateur \
    -profile slurmSingularity \
    -resume \
    --n_value $Neur \
    --plink_mem_req 50GB \
    --max_plink_cores 10 \
    --output_dir all/hiv_eur \
    -w /path/to/prs/work \
    --cut_maf 0.01 \
    --addpcs 5 \
    --type_sumstat other \
    --sumstat_head=eur \
    --sumstat_head_chr chr.38 \
    --sumstat_head_bp bp.38 \
    --sumstat_head_rs RSID \
    --sumstat_head_A1 allele1.38 \
    --sumstat_head_A2 allele0.38 \
    --sumstat_head_beta BETA.hg38 \
    --sumstat_head_se SE \
    --sumstat_head_pval P \
    --check_infoprsc 0 \
    --prscx 0 \
    --snpinfo_prcsx snpinfo_mult_1kg_hm3 \
    -c checkoption.nf \
    --prsice_binary_target T \
    --used_pval_z 0 \
    --prscs 1 \
    --covariates age,sex,Study.Site \
    --prsice_all_score 0
```

The PRS workflow produces the PRS profile used in downstream analysis.

---

## `build_pcs.r`

This script calculates principal components for population-structure adjustment.

PCs are calculated separately for:

- Each study site
- Men
- Women

The genotype data are first LD-pruned and variants with MAF below 0.01 are excluded. Twenty PCs are then calculated using PLINK2.

### Code

```r
library(pROC)

data1 <- read.table(
    'all/hiv_eur/PRScs/eur/eur.profile',
    header = T
)

data3 <- read.table(
    '../data/pheno/allsite_5_format.pheno',
    header = T
)

all <- merge(data1, data3, by = c(1,2))
all$sumscore <- all$SCORE * all$CNT

dirout <- 'pcs_res/'
dir.create(dirout)

cmt <- 1
listn <- c()

for(site in unique(all$Study.Site)){
    site_all <- all[all$Study.Site == site,]
    fileouth <- paste(dirout, site, '_all', sep='')

    write.table(
        site_all[,c(1,2)],
        file = paste(fileouth, '.ind', sep=''),
        row.names = F,
        col.names = F,
        quote = F,
        sep = '\t'
    )

    system(paste(
        '~/bin/plink2 -bfile array_plk --keep ',
        fileouth,
        '.ind --indep-pairwise 100 10 0.1 -maf 0.01 --out ',
        fileouth,
        sep=''
    ))

    system(paste(
        '~/bin/plink2 -bfile array_plk --keep ',
        fileouth,
        '.ind --extract ',
        fileouth,
        '.prune.in --pca 20 -out ',
        fileouth,
        sep=''
    ))
}

all_men <- all[all$sex == 0,]

for(site in unique(all_men$Study.Site)){
    site_all <- all_men[all_men$Study.Site == site,]
    fileouth <- paste(dirout, site, '_men', sep='')

    write.table(
        site_all[,c(1,2)],
        file = paste(fileouth, '.ind', sep=''),
        row.names = F,
        col.names = F,
        quote = F,
        sep = '\t'
    )

    system(paste(
        '~/bin/plink2 -bfile array_plk --keep ',
        fileouth,
        '.ind --indep-pairwise 100 10 0.1 -maf 0.01 --out ',
        fileouth,
        sep=''
    ))

    system(paste(
        '~/bin/plink2 -bfile array_plk --keep ',
        fileouth,
        '.ind --extract ',
        fileouth,
        '.prune.in --pca 20 -out ',
        fileouth,
        sep=''
    ))
}

all_women <- all[all$sex == 1,]

for(site in unique(all_women$Study.Site)){
    site_all <- all_women[all_women$Study.Site == site,]
    fileouth <- paste(dirout, site, '_women', sep='')

    write.table(
        site_all[,c(1,2)],
        file = paste(fileouth, '.ind', sep=''),
        row.names = F,
        col.names = F,
        quote = F,
        sep = '\t'
    )

    system(paste(
        '~/bin/plink2 -bfile array_plk --keep ',
        fileouth,
        '.ind --indep-pairwise 100 10 0.1 -maf 0.01 --out ',
        fileouth,
        sep=''
    ))

    system(paste(
        '~/bin/plink2 -bfile array_plk --keep ',
        fileouth,
        '.ind --extract ',
        fileouth,
        '.prune.in --pca 20 -out ',
        fileouth,
        sep=''
    ))
}
```

The PC files are written to:

```text
pcs_res/
```

---

## `extract_genetics_data.bash`

Bash script used to extract and prepare the PRS results generated by the PRS workflow for downstream analysis.

The PRS output is subsequently used by `extract_stat.r`.

### Code

```bash
# PRS extraction commands
# Configure the input and output paths for your local environment.

# Example:
# plink -bfile /path/to/genotype_dataset \
#     --score /path/to/prs_score_file \
#     2 4 6 \
#     --out /path/to/output
```

The exact extraction command should be configured according to the PRS output generated by the local PRS workflow.

---

## `extract_stat.r`

Main statistical analysis script.

The script:

1. Reads the PRS profile.
2. Reads phenotype data.
3. Merges PRS and phenotype information.
4. Calculates the PRS score.
5. Merges study-specific PCs.
6. Adjusts the PRS for age, sex and PCs.
7. Calculates ROC/AUC.
8. Performs logistic regression.
9. Repeats the analysis separately for men and women.
10. Writes the statistical results to CSV files.

### PRS calculation

```r
library(pROC)

data1 <- read.table(
    'all/hiv_eur/PRScs/eur/eur.profile',
    header = T
)

data2 <- read.table(
    '../data/pheno/allsite.tsv',
    header = T
)

all <- merge(
    data1,
    data2,
    by = c(1,2)
)

all$sumscore <- all$SCORE * all$CNT
```

### Number of PCs by study site

```r
pcsn = list(
    "jcs_breast" = 6,
    "jcs_oscc" = 6,
    "jcs_control" = 6,
    "jcs_cervical" = 6,
    "South" = 5,
    "ark" = 2,
    "East" = 5,
    "West" = 5
)
```

### Overall analysis

For each study site, the corresponding PC file is merged with the PRS and phenotype data:

```r
for(site in unique(all$Study.Site)){

    site_all <- all[all$Study.Site == site,]

    pcs <- read.table(
        paste('pcs_res/', site, '_all.eigenvec', sep=''),
        comment.char = '',
        header = T
    )

    site_all <- merge(
        site_all,
        pcs,
        by = c(1,2)
    )

    formula1 <- paste(
        'sumscore~',
        paste(
            c(
                'sex',
                'age',
                paste('PC', 1:pcsn[[site]], sep='')
            ),
            collapse='+'
        ),
        sep=''
    )

    site_all$sumscore_residuals <-
        glm(
            as.formula(formula1),
            data = site_all
        )$residuals

    rocauc <- roc(
        site_all$hiv,
        site_all$sumscore_residuals,
        direction = '<'
    )

    ciauc <- ci(rocauc)

    aucres <- data.frame(
        auc = rocauc$auc,
        auc_ci5 = ciauc[1],
        auc_ci95 = ciauc[2]
    )

    formula2 <- paste(
        'hiv~sumscore+',
        paste(
            c(
                'sex',
                'age',
                paste('PC', 1:pcsn[[site]], sep='')
            ),
            collapse='+'
        ),
        sep=''
    )

    gm1 <- glm(
        as.formula(formula2),
        data = site_all,
        family = "binomial"
    )

    sm <- summary(gm1)
    sm2 <- sm$coefficients[2,]

    if(cmt == 1)
        aucresf <- aucres
    else
        aucresf <- rbind(aucresf, aucres)

    listn <- c(listn, nrow(site_all))

    if(cmt == 1)
        smf <- sm2
    else
        smf <- rbind(smf, sm2)

    cmt <- cmt + 1
}
```

### Overall output

```r
smf <- cbind(
    study.site = unique(all$Study.Site),
    as.data.frame(smf),
    aucresf,
    n = listn
)

write.csv(
    smf,
    file = 'all_prsrs.csv',
    row.names = F
)
```

### Men

The male analysis uses:

```r
all_men <- all[all$sex == 0,]
```

The PRS is adjusted for age and study-specific PCs:

```r
site_all_men$sumscore_residuals <-
    glm(
        as.formula(formula1),
        data = site_all_men
    )$residuals
```

ROC/AUC and logistic regression are then performed for each study site.

The results are written to:

```text
all_men_prsrs.csv
```

### Women

The female analysis uses:

```r
all_women <- all[all$sex == 1,]
```

The PRS is adjusted for age and study-specific PCs.

The results are written to:

```text
all_women_prsrs.csv
```

---

## `checkval.r`

Script used for validation/checking of the PRS analysis.

The validation script uses the PRS and phenotype data and performs the corresponding statistical checks.

### Code

```r
library(pROC)

data1 <- read.table(
    'random1/hiv_eur/PRScs/eur/eur.profile',
    header = T
)

data2 <- read.table(
    '../data/pheno/allsite.tsv',
    header = T
)

data3 <- read.table(
    'random1/hiv_eur/pheno/allsite_5_format.pheno',
    header = T
)

all <- merge(data1, data2, by = c(1,2))
all <- merge(all, data3, by = c(1,2))

all$sumscore <- all$SCORE * all$CNT

cmt <- 1
listn <- c()

for(site in unique(all$Study.Site)){

    site_all <- all[all$Study.Site == site,]

    site_all$sumscore_residuals <-
        glm(
            sumscore ~ Pcs_1 + Pcs_2 + Pcs_3 + Pcs_4 + sex + age,
            data = site_all
        )$residuals

    rocauc <- roc(
        site_all$hiv,
        site_all$sumscore_residuals,
        direction = '<'
    )

    ciauc <- ci(rocauc)

    aucres <- data.frame(
        auc = rocauc$auc,
        auc_ci5 = ciauc[1],
        auc_ci95 = ciauc[2]
    )

    gm1 <- glm(
        hiv ~ sumscore + sex + age +
            Pcs_1 + Pcs_2 + Pcs_3 + Pcs_4,
        data = site_all
    )

    sm <- summary(gm1)
    sm2 <- sm$coefficients[2,]

    if(cmt == 1)
        aucresf <- aucres
    else
        aucresf <- rbind(aucresf, aucres)

    listn <- c(listn, nrow(site_all))

    if(cmt == 1)
        smf <- sm2
    else
        smf <- rbind(smf, sm2)

    cmt <- cmt + 1
}

smf <- cbind(
    study.site = unique(all$Study.Site),
    as.data.frame(smf),
    aucresf,
    n = listn
)

write.csv(
    smf,
    file = 'all_prsrs.csv',
    row.names = F
)
```

---

# Statistical analysis

## Logistic regression

The main association model is:

```text
HIV status ~ PRS + age + sex + principal components
```

For sex-stratified analyses, the analysis is performed within each sex group and the corresponding sex variable is not required as a covariate.

## ROC / AUC

ROC analysis evaluates the discriminatory performance of the adjusted PRS.

The analysis provides:

- AUC
- AUC confidence interval

---

# Outputs

The main statistical result files are:

```text
all_prsrs.csv
all_men_prsrs.csv
all_women_prsrs.csv
```

---

# Requirements

The workflow requires:

- R
- R package `pROC`
- PLINK
- PLINK2
- Nextflow
- PRS-CS
- Singularity/Apptainer
- A suitable computational environment for genomic analysis

Install `pROC` in R with:

```r
install.packages("pROC")
```

---

# Input data

The analysis requires:

```text
Genotype data
Phenotype data
HIV GWAS summary statistics
Principal component files
```

Individual-level genotype and phenotype data are not distributed with this repository.

Users are responsible for ensuring that they have the appropriate permissions and ethical approvals for the data they analyse.

---

# Reproducibility

The scripts are provided to document and reproduce the PRS analysis workflow.

Before running the scripts, update the input and output paths for the local environment.

In particular, the following paths should be configured:

```text
/path/to/genotype_dataset
/path/to/phenotype_data
/path/to/HIV_GWAS_summary_statistics
/path/to/PRScs_reference
```

The PRS workflow itself is maintained in:

```text
jeantristanb/qc_gwas_genet
```


---

# Citation

If you use this code or workflow, please cite the associated publication and the methods/software used for PRS calculation, genotype processing, principal component analysis, and statistical analysis.

---

# Contact

For questions, issues, or further information about the analysis, please use the GitHub repository issue tracker.
