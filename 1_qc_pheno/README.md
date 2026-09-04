# Phenotype Preparation

This repository contains the phenotype preparation for the three study datasets used in the HIV-1 GWAS analyses.

### `JCS_Pheno.txt`

JCS phenotype preparation, including:

* removal of duplicate and mistaken-identity samples;
* definition of HIV-1 status using `gag` and `env` measurements;
* addition of age and GWAS identifiers (`FID`, `IID`).

The HIV-1 positive definition is `gag >= 350` and `env >= 170`.

### `ark.md`

ARK phenotype preparation, including:

* HIV-status definition based on previous HIV testing and/or PCR;
* ancestry/PCA and relatedness QC;
* harmonisation of sex and age;
* generation of overall, female and male phenotype files.

### `awigen.md`

AWI-Gen phenotype preparation, including:

* integration of phenotype and ancestry/admixture information;
* study-site and geographic-region assignment;
* HIV-status recoding;
* ancestry and relatedness QC;
* generation of overall and sex-specific phenotype files.

Individual-level phenotype and genotype data are not included in the public repository.

