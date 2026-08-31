# HIV-1 GWAS Meta-analysis

Meta-analysis of HIV-1 GWAS summary statistics from African cohorts and an independent European-ancestry dataset.

The workflow consists of:

1. Conversion of Uganda summary statistics from GRCh37 to GRCh38
2. Preparation of summary statistics for variant matching
3. South Africa meta-analysis
4. East Africa meta-analysis
5. Sub-Saharan Africa meta-analysis
6. Sub-Saharan Africa + European-ancestry meta-analysis

Meta-analyses were performed using **METAL**.

---

# 1. Uganda summary-statistics conversion

The Uganda summary statistics were originally provided in **GRCh37** coordinates.

The conversion workflow:

* Reads compressed (`.gz`) summary-statistics files
* Converts the input format using `bcftools +munge`
* Performs GRCh37 → GRCh38 liftover using `bcftools +liftover`
* Sorts and writes the converted variants as BCF
* Formats the summary statistics using `MungeSumstats`
* Converts the summary statistics to GRCh38 coordinates
* Writes the final tab-delimited summary-statistics files

The GRCh37 and GRCh38 reference genomes and the GRCh37-to-GRCh38 chain file are required.

## 1.1 Convert compressed files to GRCh38 BCF

```bash
for File in `ls *.gz`
do
    out=${File%%.*}
    echo $out

    if [ ! -f $out".bcf" ]
    then
        exit
    fi

    zcat $File | \
    bcftools +munge \
        --no-version \
        -Ou \
        -C /dataC/alcohol/v3/Analyse/expr/bctools_test/colheaders.tsv \
        -f $HOME/GRCh37/human_g1k_v37.fasta \
        -s tmp | \
    bcftools +liftover \
        --no-version \
        -Ou \
        -- \
        -s $HOME/GRCh37/human_g1k_v37.fasta \
        -s $HOME/GRCh37/human_g1k_v37.fasta \
        -f $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
        -c $HOME/GRCh38/hg19ToHg38.over.chain.gz \
        --reject $out.reject.bcf | \
    bcftools sort \
        -o $out".bcf" \
        -Ob \
        --write-index &
done

wait
```

## 1.2 Format and liftover Uganda summary statistics

The Uganda summary statistics contain the following original fields:

```text
CHROM
GENPOS
ID
ALLELE0
ALLELE1
A1FREQ
N
TEST
BETA
SE
CHISQ
LOG10P
EXTRA
```

These are converted to a standardised format before meta-analysis.

```r
library(data.table)

options(timeout = 2000)

library(MungeSumstats)

listf <- c(
    "CHROM",
    "GENPOS",
    "ID",
    "ALLELE1",
    "ALLELE2",
    "A1FREQ",
    "N",
    "TEST",
    "BETA",
    "SE",
    "CHISQ",
    "LOG10P"
)

listf2 <- c(
    "CHR",
    "POS",
    "MarkerID",
    "A1",
    "A2",
    "EAF",
    "N",
    "TEST",
    "Beta",
    "SE",
    "CHISQ",
    "LOG10P"
)

File <- dir(pattern = "*.gz")[1]

for (File in dir(pattern = "*.gz")) {

    Data <- fread(File)

    Data2 <- Data[, ..listf]
    names(Data2) <- listf2

    Data2$Pval <- (10 ** (-Data2$LOG10P))

    MGS2 <- MungeSumstats::format_sumstats(
        Data2,
        ref_genome = "GRCh37",
        return_data = TRUE,
        allele_flip_drop = FALSE,
        bi_allelic_filter = FALSE,
        flip_frq_as_biallelic = TRUE,
        sort_coordinates = TRUE,
        imputation_ind = TRUE
    )

    MGS <- standardise_header(
        MGS2,
        return_list = FALSE
    )

    sumstats_dt_hg38 <- MungeSumstats::liftover(
        MGS,
        ref_genome = "GRCh37",
        convert_ref_genome = "GRCh38"
    )

    newfile <- paste(
        "../",
        gsub(".gz", "", File),
        sep = ""
    )

    sumstats_dt_hg38$CHR <- paste(
        "chr",
        sumstats_dt_hg38$CHR,
        sep = ""
    )

    write.table(
        sumstats_dt_hg38,
        quote = FALSE,
        row.names = FALSE,
        sep = "\t",
        file = newfile,
        col.names = TRUE
    )
}
```

The resulting Uganda summary-statistics files are then used in the meta-analysis preparation step.

---

# 2. Prepare summary statistics

A `CHRPOS` column is added to each summary-statistics file to allow variants to be matched consistently across datasets.

```bash
INPUT_DIR="summary_statistics"
OUTPUT_DIR="prepared"

mkdir -p "$OUTPUT_DIR"
```

## JCS

```bash
awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/jcs_female_breast_cancer.saige" \
> "$OUTPUT_DIR/jcs_female_breast_cancer.txt"

awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/jcs_female_cervical_cancer.saige" \
> "$OUTPUT_DIR/jcs_female_cervical_cancer.txt"

awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/jcs_female_oesophageal_cancer.saige" \
> "$OUTPUT_DIR/jcs_female_oesophageal_cancer.txt"

awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/jcs_female_non_cancer.saige" \
> "$OUTPUT_DIR/jcs_female_non_cancer.txt"

awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/jcs_male_oesophageal_cancer.saige" \
> "$OUTPUT_DIR/jcs_male_oesophageal_cancer.txt"

awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/jcs_male_non_cancer.saige" \
> "$OUTPUT_DIR/jcs_male_non_cancer.txt"
```

## ARK

```bash
awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/ark_female.saige" \
> "$OUTPUT_DIR/ark_female.txt"

awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/ark_male.saige" \
> "$OUTPUT_DIR/ark_male.txt"
```

## AWI-Gen

```bash
awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/awigen_south_female.saige" \
> "$OUTPUT_DIR/awigen_south_female.txt"

awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/awigen_south_male.saige" \
> "$OUTPUT_DIR/awigen_south_male.txt"

awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/awigen_east_female.saige" \
> "$OUTPUT_DIR/awigen_east_female.txt"

awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/awigen_east_male.saige" \
> "$OUTPUT_DIR/awigen_east_male.txt"
```

## Uganda

```bash
awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/uganda_female.txt" \
> "$OUTPUT_DIR/uganda_female_prepared.txt"

awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/uganda_male.txt" \
> "$OUTPUT_DIR/uganda_male_prepared.txt"
```

## European ancestry

```bash
awk 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,"CHRPOS",substr($0,index($0,$3))}
NR>1  {print $1,$2,$1":"$2,substr($0,index($0,$3))}' \
"$INPUT_DIR/duarte_european.txt" \
> "$OUTPUT_DIR/duarte_european_prepared.txt"
```

---

# 3. South Africa meta-analysis

The South Africa meta-analysis combines JCS, ARK and AWI-Gen South summary statistics.

```bash
cat <<EOF > meta_hiv_sa_script.txt

SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON

MARKER CHRPOS
ALLELE Allele2 Allele1
EFFECT BETA
STDERR SE
PVALUE p.value
CHROMOSOMELABEL CHR
POSITIONLABEL POS
FREQLABEL AF_Allele2
WEIGHTLABEL N_total

# JCS
PROCESS meta_fem_bc.txt
PROCESS meta_fem_cc.txt
PROCESS meta_fem_oc.txt
PROCESS meta_male_oc.txt
PROCESS meta_fem_nc.txt
PROCESS meta_male_nc.txt

# ARK
PROCESS meta_ark_fem.txt
PROCESS meta_ark_male.txt

# AWI-Gen South
PROCESS meta_awi_south_fem.txt
PROCESS meta_awi_south_male.txt

OUTFILE meta_hiv_sa .tbl
ANALYZE HETEROGENEITY
QUIT

EOF

metal meta_hiv_sa_script.txt
```

---

# 4. East Africa meta-analysis

The East Africa meta-analysis combines AWI-Gen East and Uganda summary statistics.

```bash
cat <<EOF > meta_hiv_ea_script.txt

SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON

MARKER CHRPOS
ALLELE Allele2 Allele1
EFFECT BETA
STDERR SE
PVALUE p.value
CHROMOSOMELABEL CHR
POSITIONLABEL POS
FREQLABEL AF_Allele2
WEIGHTLABEL N_total

# AWI-Gen East
PROCESS meta_awi_east_fem.txt
PROCESS meta_awi_east_male.txt

# Uganda
PROCESS meta_uganda_fem.txt
PROCESS meta_uganda_male.txt

OUTFILE meta_hiv_ea .tbl
ANALYZE HETEROGENEITY
QUIT

EOF

metal meta_hiv_ea_script.txt
```

---

# 5. Sub-Saharan Africa meta-analysis

The Sub-Saharan Africa (SSA) meta-analysis combines JCS, ARK, AWI-Gen South, AWI-Gen East and Uganda.

```bash
cat <<EOF > meta_hiv_ssa_script.txt

SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON

MARKER CHRPOS
ALLELE Allele2 Allele1
EFFECT BETA
STDERR SE
PVALUE p.value
CHROMOSOMELABEL CHR
POSITIONLABEL POS
FREQLABEL AF_Allele2
WEIGHTLABEL N_total

# JCS
PROCESS meta_fem_bc.txt
PROCESS meta_fem_cc.txt
PROCESS meta_fem_oc.txt
PROCESS meta_male_oc.txt
PROCESS meta_fem_nc.txt
PROCESS meta_male_nc.txt

# ARK
PROCESS meta_ark_fem.txt
PROCESS meta_ark_male.txt

# AWI-Gen South
PROCESS meta_awi_south_fem.txt
PROCESS meta_awi_south_male.txt

# AWI-Gen East
PROCESS meta_awi_east_fem.txt
PROCESS meta_awi_east_male.txt

# Uganda
PROCESS meta_uganda_fem.txt
PROCESS meta_uganda_male.txt

OUTFILE meta_hiv_ssa .tbl
ANALYZE HETEROGENEITY
QUIT

EOF

metal meta_hiv_ssa_script.txt
```

---

# 6. Sub-Saharan Africa + European ancestry meta-analysis

The final meta-analysis combines the African datasets with the European-ancestry summary statistics from the Duarte dataset.

```bash
cat <<EOF > meta_hiv_ssa_eur_script.txt

SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON

MARKER CHRPOS
ALLELE Allele2 Allele1
EFFECT BETA
STDERR SE
PVALUE p.value
CHROMOSOMELABEL CHR
POSITIONLABEL POS
FREQLABEL AF_Allele2
WEIGHTLABEL N_total

# JCS
PROCESS meta_fem_bc.txt
PROCESS meta_fem_cc.txt
PROCESS meta_fem_oc.txt
PROCESS meta_male_oc.txt
PROCESS meta_fem_nc.txt
PROCESS meta_male_nc.txt

# ARK
PROCESS meta_ark_fem.txt
PROCESS meta_ark_male.txt

# AWI-Gen South
PROCESS meta_awi_south_fem.txt
PROCESS meta_awi_south_male.txt

# AWI-Gen East
PROCESS meta_awi_east_fem.txt
PROCESS meta_awi_east_male.txt

# Uganda
PROCESS meta_uganda_fem.txt
PROCESS meta_uganda_male.txt

# European ancestry
PROCESS meta_duarte.txt

OUTFILE meta_hiv_ssa_eur .tbl
ANALYZE HETEROGENEITY
QUIT

EOF

metal meta_hiv_ssa_eur_script.txt
```

---

# Meta-analysis settings

All meta-analyses used the same METAL framework:

```text
Scheme:             STDERR
Marker:             CHRPOS
Effect:             BETA
Standard error:     SE
P-value:            p.value
Alleles:            Allele2 / Allele1
Frequency:          AF_Allele2
Weight:             N_total
Average frequency:  ON
Frequency range:    ON
Heterogeneity:      ANALYZE HETEROGENEITY
```

## Outputs

The main meta-analysis results are:

```text
meta_hiv_sa.tbl
meta_hiv_ea.tbl
meta_hiv_ssa.tbl
meta_hiv_ssa_eur.tbl
```

Individual-level participant data are not included in this repository. Only summary statistics that are permitted for sharing should be deposited or distributed with the analysis code.

