#!/bin/bash
#SBATCH --partition=standard
#SBATCH --time=2-00:00:00
#SBATCH --output=./popsizeABC/res/LogFiles/T1_allo_simul_data_L2_segm100_sims500_%a.log
#SBATCH -c 1 --mem=20G
#SBATCH -J T1_Allo_pABC
#SBATCH -a 1-90

cd ./popsizeABC

conda activate popsizeABC

#Original command
#python ./popsizeABC/popsizeABC_files/comp_stat1.2/simul_data.py 2 100 6250 01

#1st argument: length of chromosome segments in million basepairs
#2nd argument: number of segements 
#3rd argument: number of simulations per segment
#4th argument: batch number

for i in {1..10}; do
  Batch=$SLURM_ARRAY_TASK_ID.$i
  python ./popsizeABC_files/comp_stat1.2/simul_data.py 2 100 500 $Batch
done
