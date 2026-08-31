# Uganda HIV-1 GWAS

Genome-wide association analyses of binary HIV-1 status were performed in participants from the Uganda cohort.

Analyses were performed using **REGENIE v3.0**. The analyses used a two-step approach, with Step 1 performed using pruned genotype data and Step 2 performed using chromosome-specific imputed BGEN files.

The phenotype used for the analyses was `Y1`.

## Analyses

Several GWAS models were performed:

* HIV-1 status adjusted for principal components (PC1–PC10)
* HIV-1 status adjusted for sex
* HIV-1 status adjusted for age and sex
* HIV-1 status adjusted for sex and PC1–PC10
* HIV-1 status adjusted for age, sex and PC1–PC10
* Female-specific HIV-1 GWAS
* Female-specific HIV-1 GWAS adjusted for PC1–PC10

---

## 1. HIV-1 GWAS – PC1–PC10

Step 1 was performed using the pruned genotype dataset, adjusting for PC1–PC10.

### Step 1

```bash
regenie \
  --step 1 \
  --bed /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/ugr_new_data_38_pruned \
  --extract /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/qc/qc_pass.snplist \
  --keep /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/qc/qc_pass.id \
  --phenoFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/hiv_phenotype.txt \
  --phenoColList Y1 \
  --covarFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/cov_pca.txt \
  --covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --bsize 1000 \
  --lowmem \
  --threads 8 \
  --bt \
  --gz \
  --out /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step1/hiv_ugr_regenie_pca
```

### Step 2

The association analysis was performed chromosome by chromosome using the imputed BGEN files. The model included PC1–PC10.

```bash
for CHR in `seq 1 22`
do

BGEN_SAMPLE_PREFIX=/mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/bgen_files/chr${CHR}
BGEN_FILE=${BGEN_SAMPLE_PREFIX}.bgen
SAMPLE_FILE=${BGEN_SAMPLE_PREFIX}.sample

regenie \
  --step 2 \
  --bgen $BGEN_FILE \
  --sample $SAMPLE_FILE \
  --phenoFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/hiv_phenotype.txt \
  --phenoColList Y1 \
  --covarFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/cov_pca.txt \
  --covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --pred /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step1/hiv_ugr_regenie_pca_pred.list \
  --bsize 400 \
  --minMAC 2 \
  --threads 8 \
  --write-samples \
  --print-pheno \
  --gz \
  --bt \
  --out /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step2_hiv/hiv_pca${CHR}

done
```

---

## 2. HIV-1 GWAS – sex

This analysis adjusted for sex only.

### Step 1

```bash
regenie \
  --step 1 \
  --bed /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/ugr_new_data_38_pruned \
  --extract /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/qc/qc_pass.snplist \
  --keep /mnt/lustre01/projects/ug_gwas/users/irene/ugr_bed_file_format/qc/qc_pass.id \
  --phenoFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/hiv_phenotype.txt \
  --phenoColList Y1 \
  --covarFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/cov_pca.txt \
  --covarColList sex \
  --bsize 1000 \
  --lowmem \
  --threads 8 \
  --bt \
  --gz \
  --out /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step1/hiv_ugr_regenie_sex
```

### Step 2

```bash
for CHR in `seq 1 22`
do

BGEN_SAMPLE_PREFIX=/mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/bgen_files/chr${CHR}
BGEN_FILE=${BGEN_SAMPLE_PREFIX}.bgen
SAMPLE_FILE=${BGEN_SAMPLE_PREFIX}.sample

regenie \
  --step 2 \
  --bgen $BGEN_FILE \
  --sample $SAMPLE_FILE \
  --phenoFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/hiv_phenotype.txt \
  --phenoColList Y1 \
  --covarFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/cov_pca.txt \
  --covarColList sex \
  --pred /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step1/hiv_ugr_regenie_sex_pred.list \
  --bsize 400 \
  --minMAC 2 \
  --threads 8 \
  --write-samples \
  --print-pheno \
  --gz \
  --bt \
  --out /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step2_hiv/hiv_sex${CHR}

done
```

---

## 3. HIV-1 GWAS – age and sex

This analysis adjusted for age and sex.

### Step 1

```bash
regenie \
  --step 1 \
  --bed /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/ugr_new_data_38_pruned \
  --extract /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/qc/qc_pass.snplist \
  --keep /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/qc/qc_pass.id \
  --phenoFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/hiv_phenotype.txt \
  --phenoColList Y1 \
  --covarFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/cov_pca.txt \
  --covarColList sex,age \
  --bsize 1000 \
  --lowmem \
  --threads 8 \
  --bt \
  --gz \
  --out /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step1/hiv_ugr_regenie_age_sex
```

### Step 2

```bash
for CHR in `seq 1 22`
do

BGEN_SAMPLE_PREFIX=/mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/bgen_files/chr${CHR}
BGEN_FILE=${BGEN_SAMPLE_PREFIX}.bgen
SAMPLE_FILE=${BGEN_SAMPLE_PREFIX}.sample

regenie \
  --step 2 \
  --bgen $BGEN_FILE \
  --sample $SAMPLE_FILE \
  --phenoFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/hiv_phenotype.txt \
  --phenoColList Y1 \
  --covarFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/cov_pca.txt \
  --covarColList sex,age \
  --pred /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step1/hiv_ugr_regenie_age_sex_pred.list \
  --bsize 400 \
  --minMAC 2 \
  --threads 8 \
  --write-samples \
  --print-pheno \
  --gz \
  --bt \
  --out /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step2_hiv/hiv_age_sex${CHR}

done
```

---

## 4. HIV-1 GWAS – sex and PC1–PC10

This analysis adjusted for sex and PC1–PC10.

### Step 1

```bash
regenie \
  --step 1 \
  --bed /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/ugr_new_data_38_pruned \
  --extract /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/qc/qc_pass.snplist \
  --keep /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/qc/qc_pass.id \
  --phenoFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/hiv_phenotype.txt \
  --phenoColList Y1 \
  --covarFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/cov_pca.txt \
  --covarColList sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --bsize 1000 \
  --lowmem \
  --threads 8 \
  --bt \
  --gz \
  --out /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step1/hiv_ugr_regenie_sex_pca
```

### Step 2

```bash
for CHR in `seq 1 22`
do

BGEN_SAMPLE_PREFIX=/mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/bgen_files/chr${CHR}
BGEN_FILE=${BGEN_SAMPLE_PREFIX}.bgen
SAMPLE_FILE=${BGEN_SAMPLE_PREFIX}.sample

regenie \
  --step 2 \
  --bgen $BGEN_FILE \
  --sample $SAMPLE_FILE \
  --phenoFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/hiv_phenotype.txt \
  --phenoColList Y1 \
  --covarFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/cov_pca.txt \
  --covarColList sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --pred /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step1/hiv_ugr_regenie_sex_pca_pred.list \
  --bsize 400 \
  --minMAC 2 \
  --threads 8 \
  --write-samples \
  --print-pheno \
  --gz \
  --bt \
  --out /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step2_hiv/hiv_sex_pca${CHR}

done
```

---

## 5. HIV-1 GWAS – age, sex and PC1–PC10

This analysis adjusted for age, sex and PC1–PC10.

### Step 1

```bash
regenie \
  --step 1 \
  --bed /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/ugr_new_data_38_pruned \
  --extract /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/qc/qc_pass.snplist \
  --keep /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/qc/qc_pass.id \
  --phenoFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/hiv_phenotype.txt \
  --phenoColList Y1 \
  --covarFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/cov_pca.txt \
  --covarColList sex,age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --bsize 1000 \
  --lowmem \
  --threads 8 \
  --bt \
  --gz \
  --out /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step1/hiv_ugr_regenie_allcov
```

### Step 2

```bash
for CHR in `seq 1 22`
do

BGEN_SAMPLE_PREFIX=/mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/bgen_files/chr${CHR}
BGEN_FILE=${BGEN_SAMPLE_PREFIX}.bgen
SAMPLE_FILE=${BGEN_SAMPLE_PREFIX}.sample

regenie \
  --step 2 \
  --bgen $BGEN_FILE \
  --sample $SAMPLE_FILE \
  --phenoFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/hiv_phenotype.txt \
  --phenoColList Y1 \
  --covarFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/cov_pca.txt \
  --covarColList sex,age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --pred /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step1/hiv_ugr_regenie_allcov_pred.list \
  --bsize 400 \
  --minMAC 2 \
  --threads 8 \
  --write-samples \
  --print-pheno \
  --gz \
  --bt \
  --out /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step2_hiv/hiv_all_cov${CHR}

done
```

---

## 6. Female-specific HIV-1 GWAS

A female-specific analysis was performed using the `--sex-specific female` option in both REGENIE steps.

### Step 1

```bash
regenie \
  --step 1 \
  --bed /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/sex_ugr_new_data_38_pruned \
  --extract /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/qc/qc_pass.snplist \
  --keep /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/qc/qc_pass.id \
  --phenoFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/hiv_phenotype.txt \
  --phenoColList Y1 \
  --sex-specific female \
  --bsize 1000 \
  --lowmem \
  --threads 8 \
  --bt \
  --gz \
  --out /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step1/hiv_ugr_regenie_female
```

### Step 2

```bash
for CHR in `seq 1 22`
do

BGEN_SAMPLE_PREFIX=/mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/new_bgen_files/output/chr${CHR}
BGEN_FILE=${BGEN_SAMPLE_PREFIX}.bgen
SAMPLE_FILE=${BGEN_SAMPLE_PREFIX}.sample

regenie \
  --step 2 \
  --bgen $BGEN_FILE \
  --sample $SAMPLE_FILE \
  --phenoFile /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/hiv_phenotype.txt \
  --phenoColList Y1 \
  --sex-specific female \
  --bsize 400 \
  --minMAC 2 \
  --threads 8 \
  --write-samples \
  --print-pheno \
  --gz \
  --bt \
  --out /mnt/lustre01/projects/ug_gwas/users/irene/regenie/ugr_bed_file_format/step2_hiv/hiv_pca_female${CHR}

done
```

---

## Software

The analyses were performed using **REGENIE v3.0**. The scripts use 1 node, 8 CPU threads, 100 GB memory and a maximum runtime of 48 hours for the SLURM jobs.

The REGENIE executable used in the scripts is:

```text
/software/regenie/3.0/regenie
```

The analyses use binary-trait mode (`--bt`) and a minimum MAC threshold of 2 in Step 2.

## Data

Individual-level genotype and phenotype data are not included in this repository.

The paths shown above correspond to the original analysis environment and should be replaced with the appropriate paths when reproducing the analysis.

