```bash
# Load Nextflow
module load nextflow

# PRScs reference LD panel
refprsc_eur=<PATH_TO_PRScs_LDBLOCK_EUR>

# PLINK binary files used as the target genotype dataset
bfile=<PATH_TO_TARGET_PLINK/all_clean>

# Phenotype and summary-statistics input files
dataI=../data/pheno/allsite.tsv
sumstateur=../external_data/HIV_acquisition_meta_sample_size1_filtered_annotated.txt.hg38.bcftools

# Phenotype variable
pheno=hiv

# Effective sample size used by PRScs
Neur=594646

# Input phenotype data
data=$dataI

# Run the PRS workflow
# ~/Cancer/qc_gwas_genet corresponds to the local checkout
# of the GitHub repository jeantristanb/qc_gwas_genet
~/Cancer/qc_gwas_genet/prs.nf \
    --ref_prscs $refprsc_eur \
    --bfile $bfile \
    --data $data \
    --pheno $pheno \
    --sumstat $sumstateur \
    -profile slurmSingularity \
    -resume \
    --n_value $Neur \
    --plink_mem_req 50GB \
    --max_plink_cores 10 \
    --output_dir all/hiv_eur \
    -w <PATH_TO_NEXTFLOW_WORK_DIRECTORY> \
    --cut_maf 0.01 \
    --addpcs 5 \
    --type_sumstat other \
    --sumstat_head=eur \
    --sumstat_head_chr chr.38 \
    --sumstat_head_bp bp.38 \
    --sumstat_head_rs RSID \
    --sumstat_head_A1 allele1.38 \
    --sumstat_head_A2 allele0.38 \
    --sumstat_head_beta BETA.hg38 \
    --sumstat_head_se SE \
    --sumstat_head_pval P \
    --check_infoprsc 0 \
    --prscx 0 \
    --snpinfo_prcsx snpinfo_mult_1kg_hm3 \
    -c checkoption.nf \
    --prsice_binary_target T \
    --used_pval_z 0 \
    --prscs 1 \
    --covariates age,sex,Study.Site \
    --prsice_all_score 0
```

