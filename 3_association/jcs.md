# JCS HIV-1 GWAS

Genome-wide association analyses of binary HIV-1 status were performed in participants from the **Johannesburg Cancer Study (JCS)**.

Analyses were performed separately by sex and cancer group using **SAIGE** through the **H3ABioNet H3AGWAS pipeline**.

All analyses were adjusted for **age and the first six principal components (PC1–PC6)**.

## Analyses

The following analyses were performed:

* Female non-cancer controls
* Female breast cancer
* Female cervical cancer
* Female oesophageal squamous cell carcinoma
* Male non-cancer controls
* Male oesophageal squamous cell carcinoma

Breast and cervical cancer analyses included female participants only.

---

## Female participants

### 1. Female non-cancer controls

GWAS of binary HIV-1 status among female participants classified as non-cancer controls.

```bash
nextflow run h3abionet/h3agwas/assoc/main.nf \
 --input_dir <JCS_HG38_ARRAY_DIR> \
 --input_pat jcs \
 --data <FEMALE_NC_PHENOTYPE_FILE> \
 --pheno HIV1 \
 --print_pca 0 \
 --output_dir <OUTPUT_DIR> \
 --output fem_nc_binary \
 --sample_snps_rel 1 \
 --saige 1 \
 --covariates Age,PC1,PC2,PC3,PC4,PC5,PC6 \
 --covariates_type 0,0,0,0,0,0,0 \
 --pheno_bin 1 \
 --do_report 0 \
 --saige_otheropt_step1 ' --includeNonautoMarkersforVarRatio=TRUE ' \
 --saige_otheropt_step2 '--minInfo=0.3 --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --is_output_moreDetails=TRUE --is_imputed_data=TRUE' \
 --list_vcf <IMPUTED_VCF_LIST> \
 --vcf_field DS \
 --grm_maf 0.01 \
 --cut_maf 0.005 \
 --vcf_minmac 5 \
 -resume \
 -process.maxForks 50 \
 -profile slurmSingularity
```

### 2. Female breast cancer

GWAS of binary HIV-1 status among female participants with breast cancer.

```bash
nextflow run h3abionet/h3agwas/assoc/main.nf \
 --input_dir <JCS_HG38_ARRAY_DIR> \
 --input_pat jcs \
 --data <FEMALE_BC_PHENOTYPE_FILE> \
 --pheno HIV1 \
 --print_pca 0 \
 --output_dir <OUTPUT_DIR> \
 --output fem_bc_binary \
 --sample_snps_rel 1 \
 --saige 1 \
 --covariates Age,PC1,PC2,PC3,PC4,PC5,PC6 \
 --covariates_type 0,0,0,0,0,0,0 \
 --pheno_bin 1 \
 --do_report 0 \
 --saige_otheropt_step1 ' --includeNonautoMarkersforVarRatio=TRUE ' \
 --saige_otheropt_step2 '--minInfo=0.3 --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --is_output_moreDetails=TRUE --is_imputed_data=TRUE' \
 --list_vcf <IMPUTED_VCF_LIST> \
 --vcf_field DS \
 --grm_maf 0.01 \
 --cut_maf 0.005 \
 --vcf_minmac 5 \
 -resume \
 -process.maxForks 50 \
 -profile slurmSingularity
```

### 3. Female cervical cancer

GWAS of binary HIV-1 status among female participants with cervical cancer.

```bash
nextflow run h3abionet/h3agwas/assoc/main.nf \
 --input_dir <JCS_HG38_ARRAY_DIR> \
 --input_pat jcs \
 --data <FEMALE_CC_PHENOTYPE_FILE> \
 --pheno HIV1 \
 --print_pca 0 \
 --output_dir <OUTPUT_DIR> \
 --output fem_cc_binary \
 --sample_snps_rel 1 \
 --saige 1 \
 --covariates Age,PC1,PC2,PC3,PC4,PC5,PC6 \
 --covariates_type 0,0,0,0,0,0,0 \
 --pheno_bin 1 \
 --do_report 0 \
 --saige_otheropt_step1 ' --includeNonautoMarkersforVarRatio=TRUE ' \
 --saige_otheropt_step2 '--minInfo=0.3 --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --is_output_moreDetails=TRUE --is_imputed_data=TRUE' \
 --list_vcf <IMPUTED_VCF_LIST> \
 --vcf_field DS \
 --grm_maf 0.01 \
 --cut_maf 0.005 \
 --vcf_minmac 5 \
 -resume \
 -process.maxForks 50 \
 -profile slurmSingularity
```

### 4. Female oesophageal squamous cell carcinoma

GWAS of binary HIV-1 status among female participants with oesophageal squamous cell carcinoma.

```bash
nextflow run h3abionet/h3agwas/assoc/main.nf \
 --input_dir <JCS_HG38_ARRAY_DIR> \
 --input_pat jcs \
 --data <FEMALE_OC_PHENOTYPE_FILE> \
 --pheno HIV1 \
 --print_pca 0 \
 --output_dir <OUTPUT_DIR> \
 --output fem_oc_binary \
 --sample_snps_rel 1 \
 --saige 1 \
 --covariates Age,PC1,PC2,PC3,PC4,PC5,PC6 \
 --covariates_type 0,0,0,0,0,0,0 \
 --pheno_bin 1 \
 --do_report 0 \
 --saige_otheropt_step1 ' --includeNonautoMarkersforVarRatio=TRUE ' \
 --saige_otheropt_step2 '--minInfo=0.3 --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --is_output_moreDetails=TRUE --is_imputed_data=TRUE' \
 --list_vcf <IMPUTED_VCF_LIST> \
 --vcf_field DS \
 --grm_maf 0.01 \
 --cut_maf 0.005 \
 --vcf_minmac 5 \
 -resume \
 -process.maxForks 50 \
 -profile slurmSingularity
```

---

## Male participants

### 5. Male non-cancer controls

GWAS of binary HIV-1 status among male participants classified as non-cancer controls.

```bash
nextflow run h3abionet/h3agwas/assoc/main.nf \
 --input_dir <JCS_HG38_ARRAY_DIR> \
 --input_pat jcs \
 --data <MALE_NC_PHENOTYPE_FILE> \
 --pheno HIV1 \
 --print_pca 0 \
 --output_dir <OUTPUT_DIR> \
 --output male_nc_binary \
 --sample_snps_rel 1 \
 --saige 1 \
 --covariates Age,PC1,PC2,PC3,PC4,PC5,PC6 \
 --covariates_type 0,0,0,0,0,0,0 \
 --pheno_bin 1 \
 --do_report 0 \
 --saige_otheropt_step1 ' --includeNonautoMarkersforVarRatio=TRUE ' \
 --saige_otheropt_step2 '--minInfo=0.3 --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --is_output_moreDetails=TRUE --is_imputed_data=TRUE' \
 --list_vcf <IMPUTED_VCF_LIST> \
 --vcf_field DS \
 --grm_maf 0.01 \
 --cut_maf 0.005 \
 --vcf_minmac 5 \
 -resume \
 -process.maxForks 50 \
 -profile slurmSingularity
```

### 6. Male oesophageal squamous cell carcinoma

GWAS of binary HIV-1 status among male participants with oesophageal squamous cell carcinoma.

```bash
nextflow run h3abionet/h3agwas/assoc/main.nf \
 --input_dir <JCS_HG38_ARRAY_DIR> \
 --input_pat jcs \
 --data <MALE_OC_PHENOTYPE_FILE> \
 --pheno HIV1 \
 --print_pca 0 \
 --output_dir <OUTPUT_DIR> \
 --output male_oc_binary \
 --sample_snps_rel 1 \
 --saige 1 \
 --covariates Age,PC1,PC2,PC3,PC4,PC5,PC6 \
 --covariates_type 0,0,0,0,0,0,0 \
 --pheno_bin 1 \
 --do_report 0 \
 --saige_otheropt_step1 ' --includeNonautoMarkersforVarRatio=TRUE ' \
 --saige_otheropt_step2 '--minInfo=0.3 --is_Firth_beta=TRUE --pCutoffforFirth=0.05 --is_output_moreDetails=TRUE --is_imputed_data=TRUE' \
 --list_vcf <IMPUTED_VCF_LIST> \
 --vcf_field DS \
 --grm_maf 0.01 \
 --cut_maf 0.005 \
 --vcf_minmac 5 \
 -resume \
 -process.maxForks 50 \
 -profile slurmSingularity
```

---

## Common parameters

All six analyses used the same GWAS settings:

* **Phenotype:** `HIV1` (binary)
* **Covariates:** Age, PC1–PC6
* **GRM MAF:** 0.01
* **Association MAF:** 0.005
* **Minimum MAC:** 5
* **Minimum imputation INFO:** 0.3
* **Genotype dosage field:** `DS`
* **Firth correction:** enabled
* **H3AGWAS:** SAIGE association workflow
* **Execution:** SLURM + Singularity

Replace the following placeholders with the appropriate files and directories:

```text
<JCS_HG38_ARRAY_DIR>
<OUTPUT_DIR>
<IMPUTED_VCF_LIST>
<FEMALE_NC_PHENOTYPE_FILE>
<FEMALE_BC_PHENOTYPE_FILE>
<FEMALE_CC_PHENOTYPE_FILE>
<FEMALE_OC_PHENOTYPE_FILE>
<MALE_NC_PHENOTYPE_FILE>
<MALE_OC_PHENOTYPE_FILE>
```

Individual-level JCS genotype and phenotype data are not included in this repository.

