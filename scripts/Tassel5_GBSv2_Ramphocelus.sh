#!/bin/bash
#SBATCH -p standard
#SBATCH -N 1 
#SBATCH --output=./GBS_Ramphocelus/logFiles/Tassel.log 
#SBATCH -c 12 --mem=125G

##### for all the commands with ./tassel5/run_pipeline.pl UPDATE THE PATHS to your path to tassel5 installation

srun hostname
echo Running on host `hostname`
echo Time is `date`

#SNP calling from Genotyping by Sequencing (GBS) data, using TASSEL5 GBSv2 pipeline.
#change all the paths to the directories where you have your files!
##This pipeline requires a reference genome. I am using my own de novo assembly called in step 3! 

module load java
java -version

#Step 1 - GBSSeqToTagDBPlugin identifies tags and the taxa from the fastq files and store in the local database.
#I have the minQS default parameter (20) because I want to filter SNPs by quality 

./GBS_Ramphocelus/tassel5/run_pipeline.pl -Xms100G -Xmx120G -fork1 -GBSSeqToTagDBPlugin -i ./GBS_Ramphocelus/TEST/fastq -k ./GBS_Ramphocelus/TEST/key/barcode_key.txt -db ./GBS_Ramphocelus/TEST/output/GBS_ramphocelus.db -e ApeKI -mnQS 20 -kmerLength 64 -minKmerL 20 -mxKmerNum 1000000000 -endPlugin -runfork1 >> ./GBS_Ramphocelus/TEST/output/ramphocelus_pipeline.out

#Step2 - TagExportToFastqPlugin retrieves distinct tags stored in the database and reformats them to a FASTQ file that can be read by the Bowtie2 or BWA aligner program. 
#This output file is input to the aligner, which creates a .sam file needed for calling SNPs further down in the GBS analysis pipeline.

./GBS_Ramphocelus/tassel5/run_pipeline.pl -Xms100G -Xmx120G -fork1 -TagExportToFastqPlugin -db ./GBS_Ramphocelus/TEST/output/GBS_ramphocelus.db -o ./GBS_Ramphocelus/TEST/output/tagsForAlign.fa.gz -c 1 -endPlugin -runfork1 >> ./GBS_Ramphocelus/TEST/output/ramphocelus_pipeline.out

#Step 3 - Alignment. Build an index using the bowtie2-build command (needed for the alignment step) with the reference genome as input parameter and a named directory as an output parameter.
#Once index has been created, the alignment command can be run

module load bowtie2

bowtie2-build -f ./GBS_Ramphocelus/TEST/referenceGenome/UROC_RaFlam.v1_Scaffolded.fasta ./GBS_Ramphocelus/TEST/referenceGenome/UROC_RaFlam.v1_Scaffolded.fasta

bowtie2 -p 12 --very-sensitive -x ./GBS_Ramphocelus/TEST/referenceGenome/UROC_RaFlam.v1_Scaffolded.fasta -U ./GBS_Ramphocelus/TEST/output/tagsForAlign.fa.gz -S ./GBS_Ramphocelus/TEST/output/tagsForAlignFullvs2.sam >> ./GBS_Ramphocelus/TEST/output/ramphocelus_pipeline.out

#Step 4 - SAMToGBSdbPlugin reads a SAM file to determine the potential positions of Tags against the reference genome.

./GBS_Ramphocelus/tassel5/run_pipeline.pl -Xms100G -Xmx120G -fork1 -SAMToGBSdbPlugin -i ./GBS_Ramphocelus/TEST/output/tagsForAlignFullvs2.sam -db ./GBS_Ramphocelus/TEST/output/GBS_ramphocelus.db -endPlugin -runfork1 >> ./GBS_Ramphocelus/TEST/output/ramphocelus_pipeline.out

#Step 5 - DiscoverySNPCallerPluginV2 takes a GBSv2 database file as input and identifies SNPs from the aligned tags. Tags positioned at the same physical location are aligned against one another, SNPs are called from the aligned tags, and the SNP position and allele data are written to the database.

./GBS_Ramphocelus/tassel5/run_pipeline.pl -Xms100G -Xmx120G -fork1 -DiscoverySNPCallerPluginV2 -db ./GBS_Ramphocelus/TEST/output/GBS_ramphocelus.db -mnLCov 0.1 -mnMAF 0.01 -deleteOldData true -endPlugin -runfork1 >> ./GBS_Ramphocelus/TEST/output/ramphocelus_pipeline.out

#Step 6 - SNPQualityProfilerPlugin This plugin scores all discovered SNPs for various coverage, depth and genotypic statistics for a given set of taxa.

./GBS_Ramphocelus/tassel5/run_pipeline.pl -Xms100G -Xmx120G -fork1 -SNPQualityProfilerPlugin -db ./GBS_Ramphocelus/TEST/output/GBS_ramphocelus.db -statFile ./GBS_Ramphocelus/TEST/output/ramphocelus_SNPqual_stats.txt -deleteOldData true -endPlugin -runfork1 >> ./GBS_Ramphocelus/TEST/output/ramphocelus_pipeline.out

#Step7 - UpdateSNPPositionQualityPlugin This plugin update the data base with quality scores

./GBS_Ramphocelus/tassel5/run_pipeline.pl -Xms100G -Xmx120G -fork1 -UpdateSNPPositionQualityPlugin -db ./GBS_Ramphocelus/TEST/output/GBS_ramphocelus.db -qsFile ./GBS_Ramphocelus/TEST/output/ramphocelus_SNPqual_stats.txt -endPlugin -runfork1 >> ./GBS_Ramphocelus/TEST/output/ramphocelus_pipeline.out

#Step 8 - ProductionSNPCallerPluginV2 This plugin converts data from fastq and keyfile to genotypes, then adds these to a genotype file in VCF or HDF5 format. VCF is the default output. 
#An HDF5 file may be requested by using the suffix ".h5" on the file used in the output file parameter. Merging of samples to the same HDF5 output file may be accomplished by using the \u2013ko option described below. 
#I have the minQS default parameter (30 or 20) because I want to see how many SNPs do I get
#The barcode_ket.txt file has only the individuals I'm interested in. 

./GBS_Ramphocelus/tassel5/run_pipeline.pl -Xms100G -Xmx120G -fork1 -ProductionSNPCallerPluginV2 -db ./GBS_Ramphocelus/TEST/output/GBS_ramphocelus.db -e ApeKI -mnQS 20 -i ./GBS_Ramphocelus/TEST/fastq/ -k ./GBS_Ramphocelus/TEST/key/barcode_key.txt -kmerLength 64 -o ./GBS_Ramphocelus/TEST/output/RaFlam.v1_AllSamples.vcf -endPlugin -runfork1 >> ./GBS_Ramphocelus/TEST/output/ramphocelus_pipeline.out

#Step 9 - GetTagSequenceFromDBPlugin To extract the Tag sequences

./GBS_Ramphocelus/tassel5/run_pipeline.pl -Xms100G -Xmx120G -fork1 -GetTagSequenceFromDBPlugin -db ./GBS_Ramphocelus/TEST/output/GBS_ramphocelus.db -o ./GBS_Ramphocelus/TEST/output/Ramphocelus.tags.txt -endPlugin -runfork1 

#Step 10 - Filter the vcf file on the visual interface calling tassel with the command: ./GBS/tassel5/start_tassel.pl and save as a VCF file in the file folder or filter using vcftools 

echo Time is `date`
