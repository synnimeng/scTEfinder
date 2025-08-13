scTEfinder
(Snakemake pipeline)
This pipeline is designed to call transposable elements (TEs) from a set of fq.gz reads.
Please check the running steps in runSMK.sh;
scTE-Conda.sh provides a complete guide for setting up the runtime environment.
Description
Project directory structure and purpose


Directory	Description
configs	Configuration files required by the pipeline. At minimum, three files must be present: tools.tsv, default.tsv, and job.tsv
data	Bundled reference data needed by the pipeline, e.g. genomes, indices, barcode subsets, etc.
modules	Modular .smk scripts for the workflow. The scripts/ sub-directory contains utility scripts (Python, R)
logs	Log-file storage. We recommend redirecting logs with cmd &> logs/task.log
inputs	Example input fq.gz files. You may store your own data elsewhere—this path can be ignored.
outputs	Example output folder. You may redirect outputs elsewhere—this path can be ignored.
Configs
How to write configuration files
tools.tsv – Tool-path mapping


Column	Meaning	Notes
tool	Tool name	Must match the default file; do not change
path	Absolute path to the executable	Configure according to your system. When migrating environments, this file usually needs updating
job.tsv & default.tsv – Job-parameter mapping
Note:
default.tsv supplies default values only when a key is missing in job.tsv.
Only the parameters listed below as “Default configurable” can be set in default.tsv.


Column	Meaning	Default configurable?
Job	Job name; also used to name output files	—
DataDir	Root directory of input files (convenience field; can be empty)	✅
R1	Path to Read 1 file (barcode). If DataDir is set, becomes {DataDir}/{R1}	—
R2	Path to Read 2 file (read). If DataDir is set, becomes {DataDir}/{R2}	—
RU	Path to Read U file (UMI). If DataDir is set, becomes {DataDir}/{RU}	—
Prefix	Ignored if R1, R2, RU are provided. Common prefix of reads; pipeline then expects {DataDir}/{Prefix}_[R1/R2/RU].{Postfix} (2–3 inputs)	—
Postfix	Ignored if R1, R2, RU are provided. Common suffix/extension of reads	✅
Platform	Sequencing platform. Allowed values: 10x3-v1, 10x3-v1.r3, 10x3-v2, 10x3-v3, 10x5, 10x5-v1, 10x5-v2	✅
Tissue	Tissue type(s) for CellTypist annotation. Multiple tissues allowed, e.g. [bone_marrow, pbmc] or [bone_marrow pbmc] (comma or space delimiter). Must be valid keys in CellTypist’s tissue dictionary	✅
RefGenome	Path to reference-genome directory	✅
RefGenomeIndex	Path to the STAR genome index file	✅
CellTypistModel	Path to the CellTypist model directory	✅
OutDir	Root output directory for this job. Results are written to {OutDir}/{Job}	✅
Memory	Memory limit (passed to any tool that supports it)	✅
Thread	Number of threads/cores (passed to any tool that supports it)	✅
CondaEnv	Path to the conda environment that supports sceasy Seurat → H5AD conversion	✅
Description	Optional free-text description of the job	✅
Outputs
Output directory contents
表格
复制
Directory	Description
0.preFq	Pre-processing of FASTQs to harmonize v1 chemistry. If no special handling is required, symlinks (ln -s) are created instead of copies. Be careful not to delete source FASTQs—these are not full copies. For v1 data, an _isalign.flag file indicates whether SeqIDs are aligned
1.Map	STAR alignment results
2.Filter	QC & filtering of STAR gene-expression matrices. High-quality cell barcodes are retained. Also contains CellTypist annotation results
3.scTE	scTE output files (csv.gz)
benchmark	Per-step runtime and CPU-usage metrics. Useful for resource optimization
log	Log files (.out and .err) capturing stdout and stderr
