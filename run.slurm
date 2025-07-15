#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=1-00:00:00
#SBATCH --output=log/main_%j
##SBATCH --mail-user=user@example.com
##SBATCH --mail-type=BEGIN,END 

module load mamba
eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run main.nf \
    -profile cluster \
    -resume \
    -with-report \
    -with-timeline \
    -ansi-log false \
    -params-file params.json