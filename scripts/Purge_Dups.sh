#!/bin/bash
#SBATCH -p standard
#SBATCH --time=5-00:00:00
#SBATCH --job-name=PD_RFI301
#SBATCH --output=PurgeDupsRFI301_manual.log
#SBATCH -c 12 --mem=120G

echo Running on host `hostname`
echo Started at `date`

source ~/.bashrc
module load minimap2 

export PATH=/scratch/mcastano/programs/purge_dups/src:$PATH
PURGEDUPS=/scratch/mcastano/programs/purge_dups

PRIMARY=./RFI301_v1.asm.p.fa
ALTERNATE=./RFI301_v1.asm.a.fa
READS=./RFI301.hifi_reads.fasta.gz

cd ./PurgeDups
mkdir Primary
mkdir Alternate

#STEP 1 - Run minimap2 to align pacbio data to both assemblies generate paf files, then produce .base.cov and .stat files with the purge_dups scripts. This is to calculate the coverage cutoffs for purge dups

cd ./PurgeDups/Primary
minimap2 -x asm20 -t 12 $PRIMARY $READS | gzip -c - > RFI301_v1_asm_p.paf.gz
$PURGEDUPS/bin/pbcstat *.paf.gz
$PURGEDUPS/bin/calcuts PB.stat > cutoffs 2>calcults.log

#STEP 2 - Split the assembly and do self alignment

$PURGEDUPS/bin/split_fa $PRIMARY > RFI301_v1_asm_p.split
minimap2 -x asm5 -DP RFI301_v1_asm_p.split RFI301_v1_asm_p.split | gzip -c - > RFI301_v1_asm_p.split.self.paf.gz

#STEP 3 - Purge haplotigs and overlaps with the following command

$PURGEDUPS/bin/purge_dups -2 -T cutoffs -c PB.base.cov RFI301_v1_asm_p.split.self.paf.gz > dups.bed 2> pri_purge_dups.log

#STEP 4 - Get purged primary and haplotig sequences from draft assembly. This command only removes haplotypic duplications at the ends of contigs. If you also want to remove duplications in the middLE remove -e option (SUPER  RISKY) may delete false positive duplications.

$PURGEDUPS/bin/get_seqs -e dups.bed $PRIMARY

HAPLOTIGS=/scratch/juy3_lab/RamphocelusHiFi/HifiasmAssemblies/RFI301/PurgeDups/Primary/hap.fa

#STEP 5 - IF CUTOFFS DIDN'T WORK AUTOMATICALLY CHECK COVERAGE HISTOGRAM AND SET YOUR OWN CUTOFFS !!! - checking is a good practice so just always check
module load miniconda3
$PURGEDUPS/scripts/hist_plot.py -c cutoffs PB.stat PB.cov.png
module unload miniconda3

echo Time is `date`

#STEP 6 - Then manually change to set low, mid and high coverages based on the histogram:
#Contigs with average coverage below the low coverage threshold are set to 'JUNK' in purge_dups' BED output.
#mid coverage represents the transition between haploid and diploid coverages. Contigs with average coverage below mid coverage are tested for haplotypic duplications.
#Contigs with average coverage above high coverage are used for classifying contigs as 'REPEAT' in purge_dups' BED output.
$PURGEDUPS/bin/calcuts -l 20 -m 61 -u 130 PB.stat > cutoffs_manual
$PURGEDUPS/bin/purge_dups -2 -T cutoffs_manual -c PB.base.cov RFI301_v1_asm_p.split.self.paf.gz > dups_manual.bed 2> pri_purge_dups_manual.log
mkdir manual_coverage
cd manual_coverage
$PURGEDUPS/bin/get_seqs -e ../dups_manual.bed $PRIMARY

#STEP 7 - Merge alt.fa and alt_asm and run the steps above to get a decent (alternate assembly)

#If doing with default coverage calculation hash the "Manual coverage settings" below
#cd ./PurgeDups/Alternate
#cat $ALTERNATE $HAPLOTIGS > RFI301_v1.asm.merged_alternate.fasta
#MERGED=./PurgeDups/Alternate/RFI301_v1.asm.merged_alternate.fasta

#If doing with manual coverage settings  hash the "default coverage settings" above
mkdir manual_coverage
cd manual_coverage
HAPLOTIGS2=./PurgeDups/Primary/manual_coverage/hap.fa
cat $ALTERNATE $HAPLOTIGS2 > RFI301_v1.asm.merged_alternate.fasta 
 
MERGED=./PurgeDups/Alternate/manual_coverage/RFI301_v1.asm.merged_alternate.fasta

#SAME STEPS AS WITH THE PRIMARY ASSEMBLY
minimap2 -x asm20 -t 12 $MERGED $READS | gzip -c - > RFI301_v1.asm.merged_alternate.paf.gz
$PURGEDUPS/bin/pbcstat RFI301_v1.asm.merged_alternate.paf.gz
$PURGEDUPS/bin/calcuts PB.stat > cutoffs 2>calcults.log
$PURGEDUPS/bin/split_fa $MERGED > RFI301_v1.asm.merged_alternate.split
minimap2 -x asm5 -DP RFI301_v1.asm.merged_alternate.split RFI301_v1.asm.merged_alternate.split | gzip -c - > RFI301_v1.asm.merged_alternate.split.self.paf.gz
$PURGEDUPS/bin/purge_dups -2 -T cutoffs -c PB.base.cov RFI301_v1.asm.merged_alternate.split.self.paf.gz > dups.bed 2> alt_purge_dups.log
$PURGEDUPS/bin/get_seqs -e dups.bed $MERGED

module load miniconda3
$PURGEDUPS/scripts/hist_plot.py -c cutoffs PB.stat PB.cov.png

echo Finished at `date`
