#example blueprint file
#input setting
popid: two-epoch_fold # id of the population (no white space)
nseq: 10 # number of sequences
L:563800 # total number of observed nucleic sites, including polymorphic and monomorphic --> TO CALCULATE THIS: module load bcftools samtools/1.7 | bgzip $FA_T2 | tabix -p vcf $FA_T2.gz | bcftools query -f '%CHROM\t%POS\n' $FA_T2.gz > positions.txt | awk 'NR==1 {prev_pos=$2; prev_chr=$1; next} {if ($1 == prev_chr) {total += ($2 - prev_pos)}; prev_pos = $2; prev_chr = $1} END {print total}' positions.txt
whether_folded: true # whether the SFS is folded (true or false)
SFS: 	2839.121212121212	596.3787878787888	351.1363636363637	264.8484848484852	238.0151515151518	121.5000000000001 # snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space) --> THE NUMBERS THAT YOU IMPUT HERE NEED TO BE nseq/2 or it will give an error. eg: 12 sequences, you input 6 numbers. 
#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment this line and change this number to 2
#largest_size_of_SFS_bin_used_for_estimation: 6 # default is nseq/2 for folded SFS
pct_training: 0.67 # percentage of sites for training
nrand: 2 5 7 10 # number of random break points for each try (separated by white space)
project_dir: two-epoch_fold_FT1 # project directory
stairway_plot_dir: stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
#random_seed: 6
#output setting
mu: 4.47e-9 # assumed mutation rate per site per generation
year_per_generation: 2.6 # assumed generation time (in years)
#plot setting
plot_title: Flammigerus_T1 # title of the plot
xrange: 0.1,1000 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,0 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size
