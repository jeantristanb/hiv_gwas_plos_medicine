# AWIGEN HIV-1 GWAS

Genome-wide association analyses of binary HIV status were performed in the AWIGEN cohort using the **H3AGWAS Nextflow pipeline with SAIGE**.

The analyses were performed separately by geographic region and sex. All analyses included **age and the first five principal components (PC1–PC5)** as covariates.

## Analyses

* East African females
* East African males
* South African females
* South African males

---

## 1. East African females

GWAS of binary HIV status in East African female participants, adjusted for age and PC1–PC5.

```bash
nextflow run h3abionet/h3agwas -r ver2-dev \
  --bfile /spaces/jeantristan/recall_awigen_ark_jcs/awigen/hg38_array/awigen \
  --data /spaces/jeantristan/hiv/data/awigen/pheno/East_Female.gender \
  --pheno hiv \
  --add_pcs 5 \
  --output_dir East_Female \
  --output East_Female \
  --sample_snps_rel 1 \
  --saige 1 \
  -profile slurmSingularity \
  --vcf_list /spaces/jeantristan/hiv/data/awigen/saige/listvcf \
  -resume \
  --covariates age \
  --covariates_type '0' \
  --pheno_type '1' \
  --saige_otheropt_step2 '--minInfo=0.3 --is_Firth_beta=TRUE --pCutoffforFirth=0.05' \
  --cut_maf_rel 0.01 \
  --cut_maf 0.005 \
  --association 1 \
  --vcf_minmac 10
```

---

## 2. East African males

GWAS of binary HIV status in East African male participants, adjusted for age and PC1–PC5.

```bash
nextflow run h3abionet/h3agwas -r ver2-dev \
  --bfile /spaces/jeantristan/recall_awigen_ark_jcs/awigen/hg38_array/awigen \
  --data /spaces/jeantristan/hiv/data/awigen/pheno/East_Male.gender \
  --pheno hiv \
  --add_pcs 5 \
  --output_dir East_Male \
  --output East_Male \
  --sample_snps_rel 1 \
  --saige 1 \
  -profile slurmSingularity \
  --vcf_list /spaces/jeantristan/hiv/data/awigen/saige/listvcf \
  -resume \
  --covariates age \
  --covariates_type '0' \
  --pheno_type '1' \
  --saige_otheropt_step2 '--minInfo=0.3 --is_Firth_beta=TRUE --pCutoffforFirth=0.05' \
  --cut_maf_rel 0.01 \
  --cut_maf 0.005 \
  --association 1 \
  -w /dataH/projects/jt/work_saige_hiv \
  --vcf_minmac 15
```

---

## 3. South African males

GWAS of binary HIV status in South African male participants, adjusted for age and PC1–PC5.

```bash
nextflow run h3abionet/h3agwas -r ver2-dev \
  --bfile /spaces/jeantristan/recall_awigen_ark_jcs/awigen/hg38_array/awigen \
  --data /spaces/jeantristan/hiv/data/awigen/pheno/South_Male.gender \
  --pheno hiv \
  --add_pcs 5 \
  --output_dir South_Male \
  --output South_Male \
  --sample_snps_rel 1 \
  --saige 1 \
  -profile slurmSingularity \
  --vcf_list /home/jeantristan/spaces/hiv/data/awigen/saige/listvcf \
  -resume \
  --covariates age \
  --covariates_type '0' \
  --pheno_type '1' \
  --saige_otheropt_step2 '--minInfo=0.3 --is_Firth_beta=TRUE --pCutoffforFirth=0.05' \
  --cut_maf_rel 0.01 \
  --cut_maf 0.005 \
  --association 1 \
  -w work_saige_hiv_southmale \
  --vcf_minmac 15 \
  -c /home/jeantristan/spaces/hiv/data/awigen/saige/config.nf
```

---

## 4. South African females

GWAS of binary HIV status in South African female participants, adjusted for age and PC1–PC5.

```bash
nextflow run h3abionet/h3agwas -r ver2-dev \
  --bfile /spaces/jeantristan/recall_awigen_ark_jcs/awigen/hg38_array/awigen \
  --data /spaces/jeantristan/hiv/data/awigen/pheno/South_Female.gender \
  --pheno hiv \
  --add_pcs 5 \
  --output_dir South_Female \
  --output South_Female \
  --sample_snps_rel 1 \
  --saige 1 \
  -profile slurmSingularity \
  --vcf_list /spaces/jeantristan/hiv/data/awigen/saige/listvcf \
  -resume \
  --covariates age \
  --covariates_type '0' \
  --pheno_type '1' \
  --saige_otheropt_step2 '--minInfo=0.3 --is_Firth_beta=TRUE --pCutoffforFirth=0.05' \
  --cut_maf_rel 0.01 \
  --cut_maf 0.005 \
  --association 1 \
  -w /dataH/projects/jt/work_saige_hiv \
  --vcf_minmac 15
```

---

## Common parameters

```text
Phenotype:             hiv
Phenotype type:        Binary
Covariate:             age
Principal components:  PC1–PC5
GWAS method:           SAIGE
Imputation INFO:       ≥ 0.3
Firth correction:      Enabled
GRM MAF:               0.01
Association MAF:       0.005
Reference genome:      GRCh38/hg38
Pipeline branch:       ver2-dev
```

The minimum MAC threshold was **10 for East African females** and **15 for East African males, South African males and South African females**, following the original analysis commands.

The H3AGWAS pipeline is directly retrieved from GitHub using the `ver2-dev` branch:

```bash
nextflow run h3abionet/h3agwas -r ver2-dev
```

Individual-level AWIGEN genotype and phenotype data are not included in this repository.

