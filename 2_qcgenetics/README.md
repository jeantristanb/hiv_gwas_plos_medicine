## Genetic QC

Genetic quality control was performed before imputation and downstream HIV-1 GWAS analyses.

The repository contains the following QC documentation:

| File                                  | Description                                                                                                                                                                                                    |
| ------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `jcs_ark_awigen_before_imputation.md` | Main genotype QC before imputation, including sample/variant QC, duplicate and relatedness filtering, sex harmonisation, batch QC, variant harmonisation, merging and preparation of VCF files for imputation. |
| `ark_awigen_admiture.md`              | Ancestry and PCA QC for ARK/AWI-Gen, including LD pruning, ADMIXTURE, PCA, ancestry filtering, PCA outlier detection and relatedness filtering.                                                                |
| `JCS_ADMIXTURE.txt`                   | JCS ancestry QC using ADMIXTURE and 1000 Genomes reference populations, performed separately for males and females.                                                                                            |

### Previously published genetic QC

The genetic QC for the Uganda dataset has already been described in the corresponding published work and is therefore not reproduced in detail here. The present documentation focuses on the additional processing and QC steps used for the datasets included in this analysis.

The QC-passed datasets were subsequently used for genotype harmonisation, imputation and downstream HIV-1 GWAS/PRS analyses.

