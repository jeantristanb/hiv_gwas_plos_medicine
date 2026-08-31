# ARK HIV-1 GWAS

Genome-wide association analyses of binary HIV status were performed in the **ARK cohort** using the **H3AGWAS pipeline with SAIGE**.

Analyses were performed separately in **female and male participants**. Both analyses included **age and two principal components (Pcs_1 and Pcs_2)** as covariates.

## Analyses

* ARK females
* ARK males

---

## 1. ARK females

GWAS of binary HIV status in ARK female participants, adjusted for age, `Pcs_1` and `Pcs_2`.

```bash
module load nextflow/22.10.4

nextflow run ~/Travail/git/h3agwas/assoc \
  --input_dir /spaces/jeantristan/recall_awigen_ark_jcs/ark/hg38 \
  --input_pat ark \
  --data /home/jeantristan/Ark/consortium/hiv/pheno/phenopcs/ark_hiv_females.ind \
  --pheno hiv \
  --output_dir ark_hiv_females \
  --output ark_hiv_females \
  --sample_snps_rel 1 \
  --saige 1 \
  -profile slurmSingularity \
  --list_vcf /home/jeantristan/Ark/consortium/hiv/listvcf \
  --saige 1 \
  -resume \
  --covariates age,Pcs_1,Pcs_2 \
  --covariates_type 0,0,0 \
  --pheno_bin 1 \
  --saige_otheropt_step2 '--minInfo=0.3 --is_Firth_beta=TRUE --pCutoffforFirth=0.05' \
  --grm_maf 0.01 \
  --cut_maf 0.005 \
  --vcf_minmac 5
```

---

## 2. ARK males

GWAS of binary HIV status in ARK male participants, adjusted for age, `Pcs_1` and `Pcs_2`.

```bash
module load nextflow/22.10.4

nextflow run ~/Travail/git/h3agwas/assoc \
  --input_dir /spaces/jeantristan/recall_awigen_ark_jcs/ark/hg38 \
  --input_pat ark \
  --data /home/jeantristan/Ark/consortium/hiv/pheno/phenopcs/ark_hiv_males.ind \
  --pheno hiv \
  --output_dir ark_hiv_males \
  --output ark_hiv_males \
  --sample_snps_rel 1 \
  --saige 1 \
  -profile slurmSingularity \
  --list_vcf /home/jeantristan/Ark/consortium/hiv/listvcf \
  --saige 1 \
  -resume \
  --covariates age,Pcs_1,Pcs_2 \
  --covariates_type 0,0,0 \
  --pheno_bin 1 \
  --saige_otheropt_step2 '--minInfo=0.3 --is_Firth_beta=TRUE --pCutoffforFirth=0.05' \
  --grm_maf 0.01 \
  --cut_maf 0.005 \
  --vcf_minmac 5
```

---

## Common parameters

```text
Phenotype:             hiv
Phenotype type:        Binary
Covariates:            age, Pcs_1, Pcs_2
GWAS method:           SAIGE
Imputation INFO:       ≥ 0.3
Firth correction:      Enabled
GRM MAF:               0.01
Association MAF:       0.005
Minimum MAC:            5
Reference genome:      GRCh38/hg38
Nextflow version:      22.10.4
```

Individual-level ARK genotype and phenotype data are not included in this repository.

