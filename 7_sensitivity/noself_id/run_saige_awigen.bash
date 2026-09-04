


#phenopcs/ark_hiv_females.ind

head=ark_hiv_females
for head in ark_hiv_females ark_hiv_males 
do
output=$head
bfile=/spaces/jeantristan/recall_awigen_ark_jcs/ark/hg38/ark
data=$PWD/ark_qc_justtest/$head".tsv"
if [ ! -f $data ]
then
echo "$data not found "
exit
fi
pheno=hiv
cov=age
covtype=0

listfile=$PWD/listvcf
if [ ! -f $listfile ]
then
ls /spaces/jeantristan/recall_awigen_ark_jcs/ark/hg38_topmed_2/*.vcf.gz >  $listfile
fi
mkdir $head 
cd $head
echo """#!/bin/bash -l
#SBATCH -J $output
#SBATCH -o $output"_%j.txt"
#SBATCH -e $output"_%j.txt"
#SBATCH -t 100:20:00
#SBATCH -c 2
#SBATCH -p batch
#SBATCH --mem=20GB
module load nextflow/24.10
nextflow run ~/Travail/git/h3agwas_update/h3agwas/main.nf --bfile $bfile\
 --data $data --pheno $pheno --add_pcs 2 \
 --output_dir $head --output $head \
 --sample_snps_rel 1 --saige 1 \
  -profile slurmSingularity \
 --vcf_list $listfile  -resume --covariates $cov --covariates_type '$covtype' --pheno_type '1'  --saige_otheropt_step2 '--minInfo=0.3  --is_Firth_beta=TRUE --pCutoffforFirth=0.05' --cut_maf_rel 0.01 --cut_maf 0.005 --association 1 --vcf_minmac 15  -process.clusterOptions='--exclude=n01,n11'
""" >run.bash

cd ../
done
