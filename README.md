# HIV-1 GWAS and Meta-analysis

**CCR5 regulatory variants and novel loci associated with susceptibility to HIV-1 infection in Sub-Saharan African populations**

This repository contains the analysis code and computational workflow for the HIV-1 genome-wide association study (GWAS), meta-analysis, fine-mapping, SUSIE/colocalisation, polygenic risk score (PRS), and sensitivity analyses.

The analysis is organised into directories corresponding to the main stages of the workflow.


## Repository structure

```text
.
├── 0_inital_call/
├── 1_qc_pheno/
├── 2_qcgenetics/
├── 3_association/
├── 4_metaanalysis/
├── 5_finemapping/
├── 6_susie_coloc/
├── 7_PRS/
├── 8_sensitivity/
├── LICENSE
└── README.md
```

---

## 0. `0_inital_call`

This directory contains the initial study and documentation files used to define the analysis and prepare the contributing datasets.

---

## 1. `1_qc_pheno`

This directory contains the **phenotype quality-control and preparation** workflows.

The analyses include:

* checking and defining HIV status;
* selecting eligible participants;
* applying phenotype-related exclusion criteria;
* preparing age and sex variables;
* generating phenotype files for downstream analyses.

Phenotype preparation is performed separately for the contributing datasets where required.

---

## 2. `2_qcgenetics`

This directory contains the **genetic quality-control and ancestry workflows**.

The analyses include:

* genotype quality control;
* sample and genotype missingness filtering;
* Hardy-Weinberg equilibrium filtering;
* duplicate and relatedness detection;
* ancestry/admixture analysis;
* principal-component analysis (PCA);
* identification of ancestry and PCA outliers;
* genotype harmonisation;
* preparation of genotype data for imputation.

The resulting QC-passed datasets are used for downstream GWAS analyses.

---

## 3. `3_association`

This directory contains the **HIV-1 GWAS association analyses**.

GWAS analyses were performed using the **H3AGWAS pipeline and SAIGE** for the relevant contributing cohorts and population/sex groups.

Binary HIV status was used as the phenotype, with covariates defined for each dataset.

The resulting GWAS summary statistics are used for downstream meta-analysis.

---

## 4. `4_metaanalysis`

This directory contains the **GWAS meta-analysis workflows**.

The analyses include:

* preparation and harmonisation of GWAS summary statistics;
* variant matching across datasets;
* meta-analysis using **METAL**;
* heterogeneity analyses;
* population-specific meta-analyses;
* Sub-Saharan African meta-analysis;
* combined Sub-Saharan African and European meta-analysis.

The meta-analyses are based on GWAS summary statistics from the contributing datasets.

---

## 5. `5_finemapping`

This directory contains the **fine-mapping analyses** used to investigate genomic regions showing evidence of association.

The analyses identify variants that may account for the observed association signals and generate credible sets for associated regions.

---

## 6. `6_susie_coloc`

This directory contains **SUSIE fine-mapping and colocalisation analyses**.

These analyses are used to further investigate associated genomic regions and assess whether genetic association signals may share a common underlying signal with regulatory or other molecular traits.

---

## 7. `7_PRS`

This directory contains the **polygenic risk score (PRS) analyses**.

The workflow includes:

* preparation of genotype and phenotype data;
* construction of PRS using external GWAS summary statistics;
* adjustment for population structure and relevant covariates;
* logistic regression analyses;
* evaluation of PRS performance using ROC/AUC analyses.

---

## 8. `8_sensitivity`

This directory contains **sensitivity analyses** used to assess the robustness of the main findings.

The analyses include:

* leave-one-cohort analyses;
* analyses excluding self-reported HIV status where applicable;
* repeat GWAS/meta-analysis analyses;
* comparison of sensitivity analyses with the primary results.

These analyses assess whether the observed associations are robust to changes in cohort inclusion or phenotype definitions.

---

# Overall workflow

```text
0_inital_call
     |
     v
1_qc_pheno
     |
     v
Phenotype QC and preparation
     |
     v
2_qcgenetics
     |
     v
Genetic QC and ancestry filtering
     |
     v
Imputation
     |
     v
3_association
     |
     v
HIV-1 GWAS
     |
     v
4_metaanalysis
     |
     v
GWAS meta-analysis
     |
     v
5_finemapping
     |
     v
Fine-mapping
     |
     v
6_susie_coloc
     |
     v
SUSIE and colocalisation
     |
     v
7_PRS
     |
     v
Polygenic risk score analyses
     |
     v
8_sensitivity
     |
     v
Robustness and sensitivity analyses
```

## Data availability

Individual-level genotype and phenotype data are not included in this repository.

The repository contains analysis scripts and computational workflows that can be used to reproduce the analyses when the appropriate authorised input data are available.

## License

The repository is distributed under the license specified in:

```text
LICENSE
```

