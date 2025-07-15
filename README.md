# breseq-pipeline
Processing reads with `trimmomatic` and running `breseq` for analysis.
## Running the pipeline
1. Clone the repository.
2. Edit the parameter file e.g. `params.json` and if necessary the slurm submission script `run_pipeline.slurm`.
3. Run the pipeline on the cluster:
`sbatch run_pipeline.slurm`
4. Or locally:
`bash run_pipeline.sh`
