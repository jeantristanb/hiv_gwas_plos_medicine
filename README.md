# HIV-1 GWAS and Meta-analysis

This repository contains the analysis code used for the HIV-1 genome-wide association study (GWAS), meta-analysis, and sensitivity analyses.

The analysis is organised into separate directories corresponding to the main stages of the workflow.

## Repository structure

```text
.
├── 1_qc_pheno/
├── 2_qcgenetics/
├── 3_association/
├── 4_metaanalysis/
├── sensitivity/
├── LICENSE
└── README.md
```

---

## 1. `1_qc_pheno`

This directory contains the **phenotype quality-control and preparation** scripts.

The scripts prepare the phenotype information required for the HIV-1 GWAS, including:

* checking HIV status;
* selecting eligible participant: s;
* applying phenotype-related exclusion criteria;
* preparing age and sex variables;
* generating the phenotype files used for association analyses.

The phenotype preparation is performed separately for the contributing datasets where required.

---

## 2. `2_qcgenetics`

This directory contains the **genetic quality-control and ancestry/QC workflow**.

The analyses include:

* genotype quality control;
* sample and genotype missingness filtering;
* Hardy-Weinberg equilibrium filtering;
* relatedness and duplicate detection;
* ancestry/admixture analysis;
* principal-component analysis (PCA);
* identification of PCA outliers;
* genotype harmonisation;
* preparation of genotype data for imputation.

The output of this stage is the QC-passed genotype dataset used for downstream imputation and GWAS.

---

## 3. `3_association`

This directory contains the **HIV-1 GWAS association analyses**.

GWAS analyses were performed using the **H3AGWAS pipeline and SAIGE**.

Analyses were performed for the relevant cohorts and population/sex groups using binary HIV status as the phenotype.

Covariates and other analysis parameters are specified in the individual analysis commands.

The resulting GWAS summary statistics are used as input for the meta-analysis.

---

## 4. `4_metaanalysis`

This directory contains the **GWAS meta-analysis workflow**.

The workflow includes:

* preparation and harmonisation of GWAS summary statistics;
* conversion of summary statistics where required;
* variant matching across datasets;
* meta-analysis using **METAL**;
* heterogeneity analysis;
* regional meta-analyses;
* Sub-Saharan African meta-analysis;
* combined African and European-ancestry meta-analysis.

The meta-analysis combines summary statistics from the contributing cohorts rather than individual-level data.

---

## 5. `sensitivity`

This directory contains **sensitivity analyses** used to assess the robustness of the main findings.

These analyses include:

* leave-one-population-out meta-analysis;
* analysis excluding self-identified information where applicable;
* repeat GWAS/meta-analysis analyses;
* comparison of sensitivity-analysis results with the primary meta-analysis.

The purpose is to determine whether the main associations are robust to changes in population inclusion or phenotype/participant definitions.

---

# Overall workflow

The main analysis follows the sequence:

```text
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
GWAS summary statistics
     |
     v
4_metaanalysis
     |
     v
Meta-analysis
     |
     v
sensitivity
     |
     v
Robustness and sensitivity analyses
```

## Data availability

Individual-level genotype and phenotype data are not included in this repository.

The repository contains the analysis scripts and computational workflow required to reproduce the analyses when the appropriate authorised input data are available.

## License

The repository is distributed under the license specified in:

```text
LICENSE
```


