#!/bin/bash
#SBATCH
#SBATCH -o snakemake.out
#SBATCH -e snakemake.err
#SBATCH -p campus
#SBATCH -c 1
#SBATCH -t 7-00:00:00

bash run_snakemake.bash
