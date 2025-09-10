# ReplicateHybridZones
Pipeline and scripts associated with study on replicate hybrid zones between subspecies of the Flame-rumped Tanager (*Ramphocelus flammigerus*)

## Replicate hybrid zones reveal the progression of trait introgression through time 
Overview of the code used in Castaño et al. (2025) to characterize and compare three replicated transects of the Flame-rumped (*R. f. flammigerus*) x Lemon-rumped (*R. f. icteronotus*)Tanager hybrid zone in the Western Andes of Colombia. 

## ⚙️ Software requirements

### De novo genome assembly and annotation

- [Meryl](https://github.com/marbl/meryl)
- [GenomeScope2](https://github.com/tbenavi1/genomescope2.0)
- [Hifiasm](https://github.com/chhylp123/hifiasm)
- [Purge_Dups](https://github.com/dfguan/purge_dups)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [BUSCO](https://gitlab.com/ezlab/busco)
- [Minimap2](https://github.com/lh3/minimap2)
- [Samtools](http://www.htslib.org/)
- [Qualimap](http://qualimap.conesalab.org/)
- [BlobToolKit](https://blobtoolkit.genomehubs.org/)
- [QUAST](https://github.com/ablab/quast)
- [MitoHiFi](https://github.com/marcelauliano/MitoHiFi)
- [Medusa](https://github.com/combogenomics/medusa)
- [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler)
- [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker)
- [STAR](https://github.com/alexdobin/STAR)
- [GeMoMa](https://www.jstacs.de/index.php/GeMoMa)


### Demultiplexing, mapping, calling variants and filtering

- [TASSEL 5](https://www.maizegenetics.net/tassel)
- [VCFtools](https://vcftools.github.io/index.html)

### Population genetics analysis

- [PLINK](https://www.cog-genomics.org/plink/)
- [EEMS](https://github.com/dipetkov/eems)
- [RStudio](https://posit.co/products/open-source/rstudio/)
- [HZAR (R package)](https://cran.r-project.org/package=hzar)
- [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2/)
- [PopSizeABC](https://forge-dga.jouy.inra.fr/projects/popsizeabc/)
- [StairwayPlot2](https://github.com/xiaoming-liu/stairway-plot-v2)
- [EasySFS](https://github.com/isaacovercast/easySFS)
- [gghybrid (R package)](https://github.com/ribailey/gghybrid)
- [GEMMA](https://github.com/genetics-statistics/GEMMA)

---

## De novo genome assembly and annotation

💻 All commands below were run on a Linux-based high-performance computing (HPC) cluster using a SLURM job scheduler. Memory and CPU requirements were adjusted per software, with SLURM directives set accordingly in each submission script.

### Genome assembly 

#### 1. Count k-mers in HiFi reads to build a k-mer database with Meryl
```
cd /path/to/raw/reads
meryl count k=21 output R1_RFI301_meryldb m64086e_230626_160908.hifi_reads.fasta.gz
meryl count k=21 output R2_RFI301_meryldb m64086e_230706_190935.hifi_reads.fasta.gz
meryl count k=21 output R3_RFI301_meryldb m64190e_230621_195312.hifi_reads.fasta.gz
meryl union-sum output RFI301_meryldb R1_RFI301_meryldb R2_RFI301_meryldb R3_RFI301_meryldb
meryl histogram RFI301_meryldb > RFI301_merylhisto.txt
```
#### 2. Estimate genome size, heterozygosity, and duplication using GenomeScope2
Upload the Meryl histogram "RFI301_merylhisto.txt" to the online version of [GenomeScope2](http://genomescope.org/genomescope2.0/) to run the model and plot. 
#### 3. Assemble the genome using Hifiasm and convert gfa output to fasta 
```
hifiasm -o ./RFI301_Hifiasm_v1.asm -t 12 --primary m64086e_230626_160908.hifi_reads.fastq.gz m64086e_230706_190935.hifi_reads.fastq.gz m64190e_230621_195312.hifi_reads.fastq.gz 
awk '/^S/{print ">"$2"\n"$3}'  RFI301_Hifiasm_v1.asm.p_ctg.gfa | fold >  RFI301_v1.asm.p.fa
awk '/^S/{print ">"$2"\n"$3}'  RFI301_Hifiasm_v1.asm.a_ctg.gfa | fold >  RFI301_v1.asm.a.fa
```
#### 4. Evaluate raw assembly quality, completeness and contiguity.
 > Visualize assembly quality and contiguity using [QUAST](https://github.com/ablab/quast) which gives summary statistics and plots. 
```
quast.py RFI301_v1.asm.p.fa RFI301_v1.asm.a.fa -o Quast
```
 > Evaluate genome completeness with BUSCO and the passeriformes_odb10 database
```
busco -i RFI301_v1.asm.p.fa -o Primary -l passeriformes_odb10 -m genome -c 12 -f
busco -i RFI301_v1.asm.a.fa -o Alternate -l passeriformes_odb10 -m genome -c 12 -f
```
#### 5. Remove duplicated haplotigs and redundant sequences using PurgeDups
use the script `Purge_Dups.sh` inside the scripts folder to run purge dups on the primary and alternate assemblies.
#### 6. Run step 4 (Evaluate assembly quality and completeness) again to see how the puging affected the assembly. 
#### 7. Assemble and annotate mitochondrial genome using MitoHiFi. Follow the documentation in [MitoHiFi](https://github.com/marcelauliano/MitoHiFi) to install all dependencies.  
- [x] Step 1:
> Set variables and PATH variable
```
export PATH=/software/blast/2.10.0+:$PATH
export PATH=/software/automake/1.16:$PATH
export PATH=/software/autoconf/2.69:$PATH
export PATH=/software/jre/1.8.0_111:$PATH
export PATH=/path/to/your/MitoFinder:$PATH #Set path to your MitoFinder installation
```
> If reads are in multiple files, concatenate reads into a single fasta file 
```
zcat m64086e_230626_160908.hifi_reads.fasta.gz m64086e_230706_190935.hifi_reads.fasta.gz m64190e_230621_195312.hifi_reads.fasta.gz > RFI301.hifi_reads.fasta.gz
```
> Set Paths to reads file, MitoHifi folder (cloned from github) and output directory for all MitoHiFi analyses 
```
READS='./RFI301.hifi_reads.fasta.gz'
MITOHIFIPATH='/path/to/your/MitoHiFi' 
OUTDIR='./MitoGenome'
```
- [x] Step 2: Find and download the most closely-related mitochondrial reference genome for your species from NCBI.
> Run directly on command line.
```
cd $OUTDIR
python ${MITOHIFIPATH}/src/findMitoReference.py --species "Ramphocelus flammigerus" --outfolder ${OUTDIR} --type "mitochondrion" --min_length 14000 -n 5
```
- [x] Step 3. Based on your output from Step 2, set the name of your closely related mito genome:
```
echo "The closest related mitogenome is: NC_062466.1 (Diglossa brunneiventris)"
MITOGENOME='NC_062466.1'
```
- [x] Step 4. Run MitoHifi 
> Outputs a final_mitogenome.fasta file 
```
echo "running against ${MITOGENOME} mito using ${READS} without -p 90
python ${MITOHIFIPATH}/src/mitohifi.py -r ${READS} -f ${MITOGENOME}.fasta -g ${MITOGENOME}.gb -t 12 -d -o 2
```
#### 8. Align reads or contigs to the purged primary assembly using Minimap2 for downstream analyses and to visualize coverage patterns
> used the concatenated reads for minimap2 and then samtools to sort and index the bam file
```
minimap2 -ax asm20 -t 12 RFI301_v1.asm.p_purged.fa RFI301.hifi_reads.fasta.gz | samtools sort -o aln_hifi_reads_to_RFI301_v1.asm.p_purged.sorted.bam -T reads.tmp

samtools index aln_hifi_reads_to_RFI301_v1.asm.p_purged.sorted.bam
```
#### 9. Visualize mapping quality and coverage using Qualimap
> Allocated 120G of memory but used 115G as the max in Java memsize (to prevent java mem errors).  
```
qualimap bamqc -bam aln_hifi_reads_to_RFI301_v1.asm.p_purged.sorted.bam -c -sd -nw 400 -hm 3 -outdir Qualimap --java-mem-size=115G
```
#### 10. Perform sequence similarity searches using BLAST against reference databases for contamination screening
> This step needs a lot of memory, or it will fail so plan accordingly (>300G)
```
export BLASTDB=/software/blast/2.10.0+/blastdb

blastn -query RFI301_v1.asm.p_purged.fa -db /software/blast/2.10.0+/blastdb/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -num_threads 12 -evalue 1e-25 > RFI301_v1.asm.p_purged.tsv
```
#### 11. Visualize and evaluate genome assembly quality and contamination with BlobToolKit
> I followed the instructions to install [blobtoolkit](https://blobtoolkit.genomehubs.org/install/) (NOT BLOBTOOLS2)
- [x] Step 1: create a BlobDir (name at the end of the command)
```
blobtools create --fasta /path/to/your/purged/assembly/RFI301_v1.asm.p_purged.fa RFI301_v1
```
- [x] Step 2: Add coverage (bam file). This command need the bam files to be indexed in .csi format (done with samtools) not .bai format 
```
samtools index -c aln_hifi_reads_to_RFI301_v1.asm.p_purged.mtDNA.sorted.bam
blobtools add --cov /path/to/your/minimap/output/aln_hifi_reads_to_RFI301_v1.asm.p_purged.sorted.bam RFI301_v1
```
- [x] Step 3: Add BUSCO scores from the purged assembly
```
blobtools add --busco /path/to/your/purged/busco/output/full_table.tsv RFI301_v1
```
- [x] Step 4: Add Blast Hits 
```
blobtools add --hits /path/to/your/blast/output/RFI301_v1.asm.p_purged.tsv RFI301_v1
```
- [x] Step 5: Go to the folder with the RFI301_v1 directory and open dataset in BlobToolKit Viewer with:
```
blobtools host `pwd`
```
> This will print:
- Starting BlobToolKit API on port #### 
- Starting BlobToolKit viewer on port #### 
- Visit http://localhost:#### to use the interactive BlobToolKit Viewer. 
> Open the viewer and Download the dataset as a .txt file

- [x] Step 5: Filter out contigs that match to any order other than aves to remove contamination. I used [seqkit](https://bioinf.shenwei.me/seqkit/) but any method works.

#### 12. Use [Medusa](https://github.com/combogenomics/medusa) to Scaffold the purged and contamination screened primary contigs into chromosomes based on synteny with *T. guttata* and *D. brunneiventris* 
```
conda activate medusa_env

./medusa_env/bin/python ./medusa_env/medusa2/medusa2standalone/scripts/medusa.py -f ./path/to/folder/with/Tguttata/assembly -f ./path/to/folder/with/DigBru/assembly -i /path/to/your/purged/assembly/RFI301_v1.asm.p_purged.fasta -o Scaffolded_Tgut_Dibru
```
#### 12. Join the Mitochondrial genome assembly to the Primary assembly
```
cat /path/to/your/purged/assembly/RFI301_v1.asm.p_purged.scaffolded.fa /path/to/your/Mitogenome/final_mitogenome.fasta > /path/to/your/final/assembly/UROC_RaFlam.v1_Scaffolded.fasta.gz
```
### Genome annotation
#### 1.Build a de novo repeat database with RepeatModeler
```
# Define paths
FASTA=/path/to/your/final/assembly/UROC_RaFlam.v1_Scaffolded.fasta
OUTDIR=/path/to/your/wd/RepeatModeler/01_repeatModeler-denovo-repeat-contig-complete.lib
mkdir -p $OUTDIR && cd $OUTDIR
# Build database
BuildDatabase -name UROC_RaFlam.v1_Scaffolded -engine ncbi $FASTA &> BuildDatabase_run1.log
# Run RepeatModeler
RepeatModeler -engine ncbi -pa 12 -database UROC_RaFlam.v1_Scaffolded &>> RepeatModeler_run1.log
```
#### 2.Mask the genome with RepeatMasker
> First use the Zebra finch repeat database available with the software
```
mkdir /path/to/your/wd/RepeatMasker/Tgut
RepeatMasker -pa 12 -gff -species Taeniopygia -dir /path/to/your/wd/RepeatMasker/Tgut $FASTA
```
> Then do a second run on the already masked genome (with zebra finch in the previous step) but now using the De novo generated repeat database 
```
mkdir /path/to/your/wd/RepeatMasker/Rfla
LIB=/path/to/your/wd/RepeatModeler/01_repeatModeler-denovo-repeat-contig-complete.lib/RM_###/consensi.fa.classified
RepeatMasker -pa 12 -gff -lib $LIB -dir /path/to/your/wd/RepeatMasker/Rfla /path/to/your/wd/RepeatMasker/Tgut/UROC_RaFlam.v1_Scaffolded.fasta.masked 
```
#### 3.If available map RNAseq reads to purged final masked assembly with STAR
- [x] Step 1: Generate genome indices for the masked assembly
```
mkdir /path/to/your/STARdir/RFI301_indices
STAR --runMode genomeGenerate --genomeDir /path/to/your/STARdir/RFI301_indices --genomeFastaFiles /path/to/your/wd/RepeatMasker/Rflam/UROC_RaFlam.v1_Scaffolded.fasta.masked.masked --runThreadN 21
```
- [x] Step 2: Align RNAseq reads with STAR 
> Note: Raw RNA-seq reads have already been processed through TrimGalore.
```
STAR --runMode alignReads \
  --genomeDir /path/to/your/STARdir/RFI301_indices \
  --readFilesIn ${READS_PATH}/RFI301/RFI301_1.trimmed.fq.gz ${READS_PATH}/RFI301/RFI301_2.trimmed.fq.gz \
  --twopassMode Basic \
  --outFileNamePrefix /path/to/your/STARdir/RFI301 \
  --runThreadN 21 \
  --readFilesCommand zcat \
  --outSAMstrandField intronMotif \
  --outSAMtype BAM SortedByCoordinate
```
#### 2.Use the mapped RNAseq reads and available gene annotations for *T. guttata* for homology based genome annotation with GeMoMa
```
export PATH=/path/to/GeMoMa/env/GEMOMA/bin:$PATH
 
cd ./Annotation
 
java -Xms5G -Xmx170G -jar /pat/to/GeMoMa/git/GeMoMa-1.9/GeMoMa-1.9.jar CLI GeMoMaPipeline \
AnnotationFinalizer.r=NO o=true \
t=/path/to/your/wd/RepeatMasker/Rflam/UROC_RaFlam.v1_Scaffolded.fasta.masked.masked \
s=own i=Tgut \
a=./ncbi_dataset/data/GCF_003957565.2/genomic.gff \
g=./ncbi_dataset/data/GCF_003957565.2/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna \
r=MAPPED ERE.s=FR_FIRST_STRAND \
ERE.m=/path/to/your/STARdir/RFI301AlignedSortedByCoord.out.bam \
threads=16 \
outdir=./GeMoMa
```
---
## Demultiplexing and mapping reads, calling variants and filtering
#### 1. Read demultiplexing, mapping and variant calling. 
Use the script `Tassel5_GBSv2_Ramphocelus.sh` inside the scripts folder to run the [TASSEL 5](https://www.maizegenetics.net/tassel) GBS pipeline which will demultiplex and map reads, do a quality filter and output a vcf file (This script is commented step by step). 
#### 2. Filter variants with VCFtools 
Follow the [Speciation Genomics](https://speciationgenomics.github.io/filtering_vcfs) GitHub to generate statistics from the vcf and decide which quality, depth and missingness filters filters to apply depending on your data. 
The script `SNP_stats.sh` inside the scripts folder will calculate the statistics for the vcf file according to the tutorial.
1. Remove all the individuals with low depth (this will also remove the blank).
> Field 3 in the .idepth file is the depth per individual.  
```
awk -F'\t' '$3 <= 1.3 {print $1}' 'SNP_Stats.idepth' > LowDepth.txt
vcftools --vcf RaFlam.v1_AllSamples.vcf --remove LowDepth.txt --recode --recode-INFO-all --out RaFlam.v1_AllSamples_iDepth 
```
2. Filter for site depth (remove very low depth or very high due to paralogs or repetitive content), no indels and only biallelic SNPs - DATASET A

```
vcftools --vcf RaFlam.v1_AllSamples_iDepth.recode.vcf --remove-indels --min-alleles 2 --max-alleles 2 --minDP 4 --min-meanDP 4 --max-meanDP 30 --maxDP 30 --recode --recode-INFO-all --out RaFlam.v1_AllSamples_iDepth_sDepth_indels_Biallelic
```
3. Filter for missing data (Use this dataset for FST, Introgress, Cline analysis freq) - DATASET B
```
vcftools --vcf RaFlam.v1_AllSamples_iDepth_sDepth_indels_Biallelic.recode.vcf --max-missing 0.75 --recode --recode-INFO-all --out RaFlam.v1_AllSamples_iDepth_sDepth_indels_Biallelic_Miss0.75
```
4. Filter for Minor Allele Frequency
```
vcftools --vcf RaFlam.v1_AllSamples_iDepth_sDepth_indels_Biallelic_Miss0.75.recode.vcf --maf 0.05 --recode --recode-INFO-all --out RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_Miss0.75_maf0.05
```
5. Thin this dataset to keep 1 SNP per rad loci (LD prunning for: PCA, ADMIXTURE) - DATASET C
```
vcftools --vcf RaFlam.v1_AllSamples_iDepth_sDepth_indels_Biallelic_Miss0.75_maf0.05.recode.vcf --recode --recode-INFO-all --thin 100 --out RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_Miss0.75_maf0.05_Thinned
```
##### DEMOGRAPHIC ANALYSIS DATASETS FILTERING - DATASET D
You want to maximize the number of SNPs that you can keep while allowing almost no missing data. So you'll have to subset individuals with missing data, remove outgroups and birds you aren't going to use and then filter for site missingness. You don't normally filter for maf since rare variants can be very informative about demographic processes.
Also, demographic analysis with FSC and SNAPP requires data to be unlinked (LD prunned, while popsizeABC requires non-LD prunned data). 

6. Use the file from step 3 (Dataset B) and subset only the populations/birds you are interested in. In my case I only want individuals from the transects and from the sister species.
```
vcftools --vcf RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_Miss0.75.recode.vcf --keep Transects_Passerini.txt --recode-INFO-all --recode --out RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_MaxMissing0.75_3TrRp
```
7. Filter out individuals with a lot of missing data so that you don't get rid of so many SNPs when filtering for site missingness.
>First you need to calculate missing data per individual in your new subset file:
```
vcftools --vcf RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_MaxMissing0.75_3TrRp.recode.vcf --missing-indv --out 3TrRp
awk -F'\t' '$5 <= 0.3 {print $1}' '3TrRp.imiss' > 3TrRp_MissingData.txt
vcftools --vcf RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_MaxMissing0.75_3TrRp.recode.vcf --keep 3TrRp_MissingData.txt --recode-INFO-all --recode --out RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_MaxMissing0.75_3TrRp_iMiss0.70
```
8. Filter out all the sites that are not present in 90% of the individuals 
```
vcftools --vcf RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_MaxMissing0.75_3TrRp_iMiss0.70.recode.vcf --recode-INFO-all --max-missing 0.9 --recode --out ./popsizeABC/RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_MaxMissing0.75_3TrRp_iMiss0.70_MaxMissing0.90
```
9. Filter out the sex chromosomes because they have different demographic histories than autosomes - ZW is haploid in half of the individuals (Use this file for popsizeABC)
> Subset this file based on the pure individuals that you want for each population (FLAM/ICT SYPM - FLAM/ICT ALLO). This txt files are based on having equal representation of individuals in all transects and individuals with low missing data.
```
vcftools --vcf RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_MaxMissing0.75_3TrRp_iMiss0.70_MaxMissing0.90.recode.vcf --not-chr W --not-chr Z --recode-INFO-all --recode --out RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_MaxMissing0.75_3TrRp_iMiss0.70_MaxMissing0.90_Autosomes
``` 
10. Remove linked sites - LD prunning or thinning (Use this file for FSC and SNAPP). 
```
vcftools --vcf RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_MaxMissing0.75_3TrRp_iMiss0.70_MaxMissing0.90_Autosomes.recode.vcf --thin 100 --recode-INFO-all --recode --out RFLAM_v5_All_iDepth_sDepth4_indels_Biallelic_MaxMissing0.75_3TrRp_iMiss0.70_MaxMissing0.90_Autosomes_Thinned
```
---
## Population genetics analysis 
### 1. Convert file to PLINK and run PCA 
### 2. ADMIXTURE
### 3. EEMS
### 4. Fst - Introgress
### 5. Environmental data & isoclines
### 6. Demographic analysis 
#### PopsizeABC
 
#### Stairwayplot2 

#### FastSIMCOAL2

### 7. Geographic clines
### 8. Genomic clines
### 9. GWAS





