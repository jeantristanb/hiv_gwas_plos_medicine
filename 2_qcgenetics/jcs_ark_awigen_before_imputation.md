# Genotype QC Before Imputation

This workflow prepares the combined genotype dataset for downstream imputation and GWAS analyses.

The initial combined PLINK dataset contained **25,792 individuals** before duplicate/relatedness filtering and subsequent genotype QC.

## Workflow

1. Preparation of DNA/sample identifiers
2. Identification of duplicate or highly related individuals
3. Identification of common individuals across studies
4. Harmonisation of sex information
5. Genotype QC of the main dataset
6. Separate QC of batch 6
7. Allele-frequency comparison and variant harmonisation
8. Merging of the main dataset with batch 6
9. Conversion of the final PLINK dataset to chromosome-specific VCF files
10. Downstream imputation

---

## Input genotype data

The initial genotype dataset is provided in PLINK binary format:

```text
/dataG/cancer/erica_plink_input_v7/02_merge_bcf_plk/allind_tmp
```

The corresponding `.fam` file contains:

```text
25,792 individuals
```

This represents the starting dataset before QC and filtering.

---

# 1. Prepare DNA/sample identifiers

The first step prepares the DNA/sample identifiers required for matching individuals across datasets.

```bash
Rscript ./extract_dnaid.r
```

---

# 2. Identify duplicate or highly related individuals

Using the raw genotype data, variants are first LD-pruned and then used to calculate pairwise relatedness.

The initial filtering uses:

* MAF ≥ 0.01
* Genotype missingness ≤ 0.01
* LD pruning with a window of 100 variants
* Step size of 20 variants
* r² threshold of 0.1

Only individuals listed in `dataI/ind.keep` are included.

```bash
bfile=/dataG/cancer/erica_plink_input_v7/02_merge_bcf_plk/allind_tmp
dirdata=allind

mkdir -p $dirdata

if [ ! -f $dirdata/duplicate.tmp ]
then

~/bin/plink \
  --autosome \
  -bfile $bfile \
  --threads 20 \
  -out $dirdata/hihrel \
  --maf 0.01 \
  --indep-pairwise 100 20 0.1 \
  --geno 0.01 \
  --keep dataI/ind.keep

~/bin/plink \
  -bfile $bfile \
  --threads 20 \
  -extract $dirdata/hihrel.prune.in \
  -out $dirdata/hihrel \
  --genome \
  --keep dataI/ind.keep

awk '{if($10>0.8)print $0}' \
  $dirdata/hihrel.genome \
  > $dirdata/duplicate.tmp

fi
```

Pairs with `PI_HAT > 0.8` are written to:

```text
allind/duplicate.tmp
```

---

# 3. Identify common individuals across studies

Individuals shared between the relevant study datasets are identified using:

```bash
Rscript ./extract_commonarkagincourt.r
```

This step is used to identify the individuals that can be matched across the study datasets.

---

# 4. Update sex information and select variants

Sex information is updated using:

```text
allind/udpate.sex
```

The genotype dataset is also restricted to the variants specified in:

```text
../cmpallfreq/withoutb6/postokeep.bed
```

```bash
~/bin/plink \
  -bfile $bfile \
  --update-sex allind/udpate.sex \
  --make-bed \
  -threads 20 \
  -out $dirdata/allind_tmp_b6_withsex \
  --keep dataI/ind.keep \
  --extract range ../cmpallfreq/withoutb6/postokeep.bed
```

---

# 5. Genotype QC of the main dataset

The main dataset is processed using the H3AGWAS QC workflow.

The QC thresholds are:

* Genotype missingness: 0.05
* Individual missingness: 0.05
* Hardy-Weinberg equilibrium P-value: 1 × 10⁻⁵

The duplicate list generated above is incorporated using `--qc_dup 1`.

Sex is used as the phenotype/case-control column and batch information is retained using `Batch`.

```bash
/home/jeantristan/Travail/git/h3agwas_update/h3agwas/main.nf \
  --qc_dup 1 \
  --bfile $dirdata/allind_tmp_b6_withsex \
  --data dataI/all.pheno \
  --out all \
  -profile slurmSingularity \
  -resume \
  --output_dir all_merge \
  --qc 1 \
  --case_control_col sex \
  --cut_geno 0.05 \
  --cut_mind 0.05 \
  --cut_hwe 0.00001 \
  --batch_col Batch
```

The QC-passed dataset is generated under:

```text
all_merge/qc/
```

with the main cleaned PLINK dataset:

```text
all_merge/qc/out_cleany
```

---

# 6. QC of batch 6

Batch 6 is processed separately before being merged with the main dataset.

The batch 6 genotype dataset is:

```bash
bfile6=/dataG/cancer/erica_plink_input_v7/02_merge_bcf_plk/process_b06/tmpplink/09_cancer1_batch06
```

## 6.1 Update batch 6 sample identifiers

The `.fam` identifiers are converted into the required format:

```bash
awk '{print $1"\t"$2"\t"$2"\t"$2}' \
  $bfile6".fam" \
  > dataI/update_b6.ind
```

---

## 6.2 Select variants and filter genotype missingness

Batch 6 is restricted to the variants in:

```text
../cmpallfreq/withoutb6/postokeep.bed
```

A genotype missingness threshold of 0.30 is applied at this stage.

```bash
~/bin/plink \
  -bfile $bfile6 \
  --extract range ../cmpallfreq/withoutb6/postokeep.bed \
  --update-ids dataI/update_b6.ind \
  --geno 0.3 \
  --out dataI/update_b6_clean \
  --make-bed
```

---

## 6.3 Update batch 6 sex information

Sex information is updated using:

```text
dataI/pheno_b6.sex
```

```bash
~/bin/plink \
  -bfile dataI/update_b6_clean \
  --make-bed \
  --update-sex dataI/pheno_b6.sex \
  --out dataI/update_b6_clean_sex
```

---

## 6.4 Genotype QC of batch 6

The H3AGWAS QC workflow is applied to batch 6.

QC thresholds:

* Genotype missingness: 0.05
* Individual missingness: 0.05
* HWE P-value: 1 × 10⁻⁵

Duplicate removal is disabled for this step using `--qc_dup 0`.

```bash
/home/jeantristan/Travail/git/h3agwas_update/h3agwas/main.nf \
  --qc_dup 0 \
  --bfile dataI/update_b6_clean_sex \
  --data dataI/pheno_b6.pheno \
  --out b6 \
  -profile slurmSingularity \
  -resume \
  --output_dir b6 \
  --qc 1 \
  --case_control_col sex \
  --cut_geno 0.05 \
  --cut_mind 0.05 \
  --cut_hwe 0.00001 \
  --batch_col Batch \
  --pheno sex
```

The QC-passed batch 6 dataset is generated under:

```text
b6/qc/
```

with:

```text
b6/qc/out_cleany
```

---

# 7. Calculate allele frequencies

Allele frequencies are calculated for both QC-passed datasets.

## Main dataset

```bash
~/bin/plink2 \
  -bfile all_merge/qc/out_cleany \
  --freq cols=chrom,pos,ref,alt,reffreq,nobs \
  -out all_merge
```

## Batch 6

```bash
~/bin/plink2 \
  -bfile b6/qc/out_cleany \
  --freq cols=chrom,pos,ref,alt,reffreq,nobs \
  -out b6
```

These frequency files are used for comparison and variant harmonisation between the datasets.

---

# 8. Harmonise batch 6 variant names

Batch 6 variant identifiers are updated using:

```text
update_y.update
```

```bash
~/bin/plink \
  -bfile b6/qc/out_cleany \
  --extract update_y.update \
  --update-name update_y.update 2 1 \
  -out b6/qc/out_cleany_freq \
  -make-bed
```

The resulting dataset is used for merging with the main dataset.

---

# 9. Merge the main dataset and batch 6

A directory is created for the final merged genotype dataset:

```bash
mkdir -p all_withb6/
```

The two QC-passed datasets are then merged:

```bash
~/bin/plink \
  -bfile all_merge/qc/out_cleany \
  --bmerge b6/qc/out_cleany_freq \
  --merge-mode 1 \
  -out all_withb6/all
```

The final merged PLINK dataset is:

```text
all_withb6/all
```

---

# 10. Convert the final genotype dataset to VCF

The final merged PLINK dataset is converted to VCF using H3AGWAS.

Chromosome-parallel conversion is enabled using:

```text
--convertinvcf_parralchro 1
```

```bash
/home/jeantristan/Travail/git/h3agwas_update/h3agwas/main.nf \
  --bfile all_withb6/all \
  --output_dir all_withb6/ \
  --output all_b6 \
  --convertinvcf 1 \
  -profile slurmSingularity \
  -resume \
  --convertinvcf_parralchro 1 \
  --file_ref_gzip /home/jeantristan/Data/human_9606/All_20180423.vcf.gz \
  --fasta ~/Ark/ark/array/calling_v0224/egt_2019/qc/qc_batch12_cleanpos/utils/data/hs37d5_chr.fa_clean.fa.gz
```

The resulting chromosome-specific VCF files are used for the subsequent imputation step.

---

# QC parameters

| Step                | Parameter              | Threshold |
| ------------------- | ---------------------- | --------: |
| Relatedness         | MAF                    |      0.01 |
| Relatedness         | Genotype missingness   |      0.01 |
| LD pruning          | Window                 |       100 |
| LD pruning          | Step                   |        20 |
| LD pruning          | r²                     |       0.1 |
| Duplicate screening | PI_HAT                 |     > 0.8 |
| Main QC             | Genotype missingness   |      0.05 |
| Main QC             | Individual missingness |      0.05 |
| Main QC             | HWE                    |  1 × 10⁻⁵ |
| Batch 6 initial QC  | Genotype missingness   |      0.30 |
| Batch 6 QC          | Genotype missingness   |      0.05 |
| Batch 6 QC          | Individual missingness |      0.05 |
| Batch 6 QC          | HWE                    |  1 × 10⁻⁵ |

# Workflow summary

```text
25,792 individuals
        |
        v
Raw PLINK genotype data
        |
        v
DNA/sample ID preparation
        |
        v
LD pruning
        |
        v
Relatedness / duplicate screening
        |
        v
Common-individual identification
        |
        v
Sex harmonisation
        |
        +-------------------------+
        |                         |
        v                         v
   Main dataset               Batch 6
        |                         |
        v                         v
       QC                        QC
        |                         |
        +-----------+-------------+
                    |
                    v
       Allele-frequency comparison
                    |
                    v
       Variant-name harmonisation
                    |
                    v
             Dataset merging
                    |
                    v
          Final PLINK dataset
                    |
                    v
       Chromosome-specific VCFs
                    |
                    v
               Imputation
```

## Data availability

The individual-level genotype and phenotype data used in this workflow are not included in the public repository.

The commands document the computational workflow used to prepare the genotype data for downstream imputation and GWAS analyses.

