# ARK -AWI-Gen Admixture and PCA QC

This step performs ancestry and principal-component (PCA) quality control before the AWI-Gen HIV-1 GWAS.

The workflow uses the `qc_gwas_genet` Nextflow pipeline and combines:

* LD pruning
* Reference-population admixture analysis
* Principal-component analysis
* PCA outlier detection
* Relatedness/PI_HAT filtering
* Ancestry proportion thresholds

The analysis uses **African ancestry (AFR)** as the reference population.

## Admixture and PCA QC

The following command runs the ancestry and PCA QC workflow:

```bash
~/nextflow /home/jeantristan/Cancer/qc_gwas_genet/main.nf \
  --bfile $bfile \
  -profile slurmSingularity \
  -resume \
  --bfile_ref $bfileref \
  --plink_indep_pairwise "500 100 0.1" \
  --plink_thin_count $nbsnps \
  --dataref $dataref \
  --headphenoref_pop Ancestrality \
  --bestk_admixture_popi 6 \
  --bestk_admixture_popref 4 \
  --data $data \
  --pop_ref AFR \
  --pheno $PhenoS \
  --cut_maf 0.005 \
  --headpheno_pop $headphenopop \
  --max_plink_core 20 \
  --maxadm_popref 0.7 \
  --maxadm_popi 0.9 \
  --max_plink_cores 20 \
  --max_admixture_cores 20 \
  --output_dir $Region \
  --output $Region \
  --workDir $dirwork/ \
  --eig_outliersigmathresh 5 \
  --max_eig_outliernbpcs 5 \
  --pi_hat_data 0.9
```

## Main parameters

### LD pruning

Variants are LD-pruned using:

```text
Window:       500
Step:         100
r² threshold: 0.1
```

```text
--plink_indep_pairwise "500 100 0.1"
```

A specified number of variants is subsequently retained using:

```text
--plink_thin_count $nbsnps
```

### Variant filtering

Variants are filtered using a minimum allele frequency of:

```text
MAF ≥ 0.005
```

```text
--cut_maf 0.005
```

### Admixture analysis

The admixture analysis uses:

```text
Reference population: AFR
Reference K:          4
Study population K:   6
```

```text
--pop_ref AFR
--bestk_admixture_popref 4
--bestk_admixture_popi 6
```

The ancestry proportion thresholds are:

```text
Reference population maximum threshold: 0.7
Study population maximum threshold:     0.9
```

```text
--maxadm_popref 0.7
--maxadm_popi 0.9
```

### PCA outlier filtering

PCA outliers are identified using a significance threshold of:

```text
5 standard deviations
```

with a maximum of:

```text
5 principal components
```

```text
--eig_outliersigmathresh 5
--max_eig_outliernbpcs 5
```

### Relatedness filtering

The PI_HAT threshold is:

```text
PI_HAT = 0.9
```

```text
--pi_hat_data 0.9
```

### Computational resources

The workflow allows up to:

```text
PLINK:      20 cores
ADMIXTURE:  20 cores
```

```text
--max_plink_cores 20
--max_admixture_cores 20
```

The pipeline is executed using the SLURM/Singularity profile:

```text
-profile slurmSingularity
```

and Nextflow resume functionality:

```text
-resume
```

## Input variables

The command uses the following user-defined variables:

```text
$bfile       Input PLINK genotype dataset
$bfileref    Reference PLINK genotype dataset
$dataref     Reference population metadata
$data        Study/sample metadata
$PhenoS      Phenotype information
$headphenopop Phenotype column for the study population
$Region      Output region/name
$nbsnps      Number of variants retained after thinning
$dirwork     Nextflow working directory
```

## QC summary

```text
Input genotype data
        |
        v
      MAF QC
     MAF ≥ 0.005
        |
        v
      LD pruning
   500 / 100 / r²=0.1
        |
        v
   Variant thinning
        |
        v
     PCA analysis
        |
        v
 PCA outlier filtering
   5 SD / up to 5 PCs
        |
        v
  Admixture analysis
        |
        v
 AFR reference population
        |
        v
 Ancestry proportion QC
        |
        v
 Relatedness / PI_HAT QC
        |
        v
Final AWI-Gen QC dataset
```

The resulting ancestry and PCA QC information is subsequently used in the AWI-Gen phenotype preparation and HIV-1 GWAS analyses.

