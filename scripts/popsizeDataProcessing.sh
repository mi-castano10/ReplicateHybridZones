#!/bin/bash

#Script to prepare data for input for popsizeABC script stat_from_vcf.py 

# To run this script you need the popfiles.txt and individual.txt files already in the popfiles folder 
# Asing each individual to any of these populations Pure_Flam_Allopatric Pure_Flam_Sympatric Pure_Ict_Allopatric Pure_Ict_Sympatric 

module load vcftools bcftools samtools/1.7 plink/1.9 R

#1. Filter the main VCF file with the ict or allo population files for each transect 
## CHANGE HERE THE TRANSECT NUMBER TO RUN IT FOR A DIFFERENT TRANSECT

MAINVCF=RaFlam.v1_AllSamples_iDepth_sDepth_indels_Biallelic_Miss0.75_OnlyFlam_mac1_3T_iMiss0.70_MaxMissing0.90_Autosomes.recode.vcf

cd ./popsizeABC/transect1/popfiles

for i in *c.txt; do for f in `ls -1 "$i" | sed 's/.txt//'`; do vcftools --vcf $MAINVCF --recode --recode-INFO-all --keep $i --out ./popsizeABC/transect1/vcfs/$f; done; done
 
#2. Instead of running the updateIDs script in the terminal we're going to implement the script in this pipeline

for i in *popfile.txt; do for f in `ls -1 "$i" | sed 's/_popfile.txt//'`; do awk 'BEGIN{FS=OFS="\t"}{print $1 OFS $1 OFS $2 OFS $1}' "$i" > ./popsizeABC/transect1/updateIDs/${f}_updateIDs_file.txt; done; done

#3. Instead of running the make_bed_plink script in the terminal we're going to include it in this pipeline. This script creates a bed, bim, fam, nosex file for each vcf. 

cd ./popsizeABC/transect1/vcfs

for vcf in *.vcf; do
    f=$(echo "$vcf" | sed 's/.recode.vcf//')
    plink --vcf "$vcf" \
          --allow-extra-chr \
          --chr-set 80 \
          --double-id \
          --make-bed \
          --update-ids ./popsizeABC/transect1/updateIDs/${f}_updateIDs_file.txt \
          --out "./popsizeABC/transect1/make_bed_plink/$f"
done

#4. Make folders for the allopatric and sympatric populations, copy the vcfs to them and then bgzip and tabix each file. 
#First ALLOPATRIC
cd ./popsizeABC/transect1/ALLO_scaffold

cp ./popsizeABC/transect1/vcfs/Pure_Flam_Allopatric.recode.vcf ./popsizeABC/transect1/ALLO_scaffold

cp ./popsizeABC/transect1/vcfs/Pure_Ict_Allopatric.recode.vcf ./popsizeABC/transect1/ALLO_scaffold

for f in ./*.vcf; do \
   bgzip -c $f > $f.gz; \
done

for f in ./*.vcf.gz; do \
   tabix $f; \
done

#5. split by scaffold for each vcf
for i in *.vcf; do for f in `ls -1 "$i" | sed 's/.recode.vcf//'`; do mkdir ${f}_vcf; done; done

for i in *.vcf; do for f in `ls -1 "$i" | sed 's/.recode.vcf//'`; do \
Rscript ./popsizeABC/popsizeABC_files/comp_stat1.2/splitVCFbyScaffold_unzip.R $i ${f}_vcf; done; done 

#6. Enter the folder for each population, then bgzip all the vcf files and then separate uncompressed .vcf files into another folder with the same name but ending in _vcf

#First bgzip all files
for i in *.vcf; do for f in `ls -1 "$i" | sed 's/.recode.vcf//'`; do mkdir $f; \
for file in ${f}_vcf/*.vcf; do bgzip -c $file > $file.gz; done; done; done

#Then transfer all the gz files to the population name folder!
for i in *.vcf; do for f in `ls -1 "$i" | sed 's/.recode.vcf//'`; do \
for gzvcf in ${f}_vcf/*.vcf.gz; do mv $gzvcf ./$f/; done; done; done

#THEN SYMPATRIC
cd ./popsizeABC/transect1/SYMP_scaffold

cp ./popsizeABC/transect1/vcfs/Pure_Flam_Sympatric.recode.vcf ./popsizeABC/transect1/SYMP_scaffold

cp ./popsizeABC/transect1/vcfs/Pure_Ict_Sympatric.recode.vcf ./popsizeABC/transect1/SYMP_scaffold

for f in ./*.vcf; do \
   bgzip -c $f > $f.gz; \
done

for f in ./*.vcf.gz; do \
   tabix $f; \
done

#5. split by scaffold for each vcf
for i in *.vcf; do for f in `ls -1 "$i" | sed 's/.recode.vcf//'`; do mkdir ${f}_vcf; done; done

for i in *.vcf; do for f in `ls -1 "$i" | sed 's/.recode.vcf//'`; do \
Rscript ./popsizeABC/popsizeABC_files/comp_stat1.2/splitVCFbyScaffold_unzip.R $i ${f}_vcf; done; done 

#6. Enter the folder for each population, then bgzip all the vcf files and then separate uncompressed .vcf files into another folder with the same name but ending in _vcf

#First bgzip all files
for i in *.vcf; do for f in `ls -1 "$i" | sed 's/.recode.vcf//'`; do mkdir $f;  \
for file in ${f}_vcf/*.vcf; do bgzip -c $file > $file.gz; done; done; done

#Then transfer all the gz files to the population name folder!
for i in *.vcf; do for f in `ls -1 "$i" | sed 's/.recode.vcf//'`; do \
for gzvcf in ${f}_vcf/*.vcf.gz; do mv $gzvcf ./$f/; done; done; done
