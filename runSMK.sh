
# app.conda
# conda activate te

cd /path/to/scTEfinder

TE_ENV='/path/to/envs/scte'
export LD_LIBRARY_PATH="${TE_ENV}/lib:$LD_LIBRARY_PATH"
export PATH="${TE_ENV}/bin:$PATH"

SMK='/path/to/snakemake'

# Display the pipeline without actually running anything
#$SMK --config job=configs/test.tsv -j 3 --dryrun

# Run the pipeline  
# -j 2 means allow at most 2 concurrent jobs; adjust according to resource limits
$SMK --config job=configs/demo.tsv -j 2 --unlock
$SMK --config job=configs/demo.tsv -j 2

# List all available pipeline target rules
# $SMK --config job=configs/demo.tsv --list-target-rules

#<<EOF
# Run only the "RunAnnoOnly" rule, forcing refresh of existing outputs
#$SMK --config job=configs/demo.tsv -j 3 --dryrun RunAnnoOnly

# Snakemake considers only file timestamps; code changes do not trigger re-run
#$SMK --config job=configs/demo.tsv -j 3 --dryrun --rerun-triggers mtime

# Manually change a file's timestamp In Linux
#stat file
#touch -t 202503010113.49 file
#EOF
