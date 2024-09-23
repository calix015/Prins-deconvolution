#!/bin/bash -l        
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --mem=490g
#SBATCH --tmp=100g
#SBATCH --time=8:00:00
#SBATCH -A prin0088
#SBATCH --mail-user=calix015@umn.edu
#SBATCH --mail-type=ALL
#SBATCH -e SM_MAST%j.err
#SBATCH -o SM_MAST%j.out
#SBATCH -p ag2tb

module load R/4.1.0

cd /home/riss/calix015/2023/prins/code

Rscript /home/riss/calix015/2023/prins/code/sigMatrix.MAST.r
