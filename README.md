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
Upload the Meryl histogram "RFI301_merylhisto.txt" to the online version of [GenomeScope2](https://github.com/tbenavi1/genomescope2.0) to run the model and plot. 
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
```
echo "running against ${MITOGENOME} mito using ${READS} without -p 90
python ${MITOHIFIPATH}/src/mitohifi.py -r ${READS} -f ${MITOGENOME}.fasta -g ${MITOGENOME}.gb -t 12 -d -o 2
```
#### 8. Align reads or contigs to the purged primary assembly using Minimap2 for downstream analyses and to visualize coverage patterns 
```
```
#### 9. Visualize mapping quality and coverage using Qualimap 
```
```
#### 10. Perform sequence similarity searches using BLAST against reference databases
> This step needs a lot of memory, or it will fail so plan accordingly (>300G)
```
```
#### 11. Visualize and evaluate genome assembly quality and contamination with BlobToolKit
> download the dataset and filter out contigs that match to any order other than aves to remove contamination
```
```
#### 12. Join the Mitochondrial genome assembly to the Primary assembly

### Genome annotation

#### 1.If available map RNAseq reads to purged final assembly with STAR
```
```
#### 2.Used the mapped reads and available annotations for closely related species for homology based genome annotation with GeMoMa
```
```
---
## Demultiplexing and mapping reads, calling variants and filtering
#### 1. Read demultiplexing, mapping and variant calling. 
Use the script `TasselV5.sh` inside the scripts folder to run the [TASSEL 5](https://www.maizegenetics.net/tassel) GBS pipeline which will demultiplex and map reads, do a quality filter and output a vcf file (This script is commented step by step). 
#### 2. Filter variants with VCFtools 
Follow the [Speciation Genomics](https://speciationgenomics.github.io/filtering_vcfs) GitHub to generate statistics from the vcf and decide which quality, depth and missingness filters filters to apply depending on your data. 
Then do a basic quality/depth filter using VCFtools  
```
vcftools --gzvcf 
```
Then filter your general dataset into 4 datasets according to the assumptions of the different analysis.
- [x] Dataset A ()
```
vcftools --gzvcf 
```
- [x] Dataset B ()
```
vcftools --gzvcf 
```
- [x] Dataset C ()
```
vcftools --gzvcf 
```
- [x] Dataset D ()
```
vcftools --gzvcf ()
```

---
## Population genetics analysis 
#### 1. 
#### 1. Convert file to PLINK and run PCA 
#### 2. ADMIXTURE
#### 3. EEMS
#### 4. Demographic analysis 
#### 5. Geographic clines
#### 6. Genomic clines
#### 7. GWAS





