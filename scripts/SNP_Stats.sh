#!/bin/bash
#SBATCH -p standard -t 08:00:00
#SBATCH -N 1 
#SBATCH --output=./GBS_Ramphocelus/TEST/logFiles/Tassel_SNPstats.log 
#SBATCH --mail-type=ALL --mail-user=mcastano@ur.rochester.edu
#SBATCH -c 4 --mem=40G

module load vcftools bcftools

# change for --gzvcf in the commands below if your vcf is compressed

VCF_IN=./GBS_Ramphocelus/TEST/output/RaFlam.v1_AllSamples.vcf.gz

mkdir /GBS_Ramphocelus/TEST/output/SNP_Stats

OUT=./GBS_Ramphocelus/TEST/output/SNP_Stats

vcftools --vcf $VCF_IN --freq2 --out $OUT --max-alleles 2

vcftools --vcf $VCF_IN --depth --out $OUT

vcftools --vcf $VCF_IN --site-mean-depth --out $OUT

vcftools --vcf $VCF_IN --missing-indv --out $OUT

vcftools --vcf $VCF_IN --missing-site --out $OUT

vcftools --vcf $VCF_IN --het --out $OUT

# Visualize the output of this data in R
#filter your vcf files according to the quality of the data that you see 
