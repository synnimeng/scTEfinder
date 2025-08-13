# scTEfinder 
__(snakemake pipeline)__

This pipeline is designed to call TEs from a set of fq.gz reads.

Please check the runSMK.sh script for execution steps
scTE-Conda.sh provides complete steps for setting up the runtime environment

## Description

Directory function description of this project

| dir | description |
| --- | --- |
| __configs__ | Configuration files required by the pipeline, __at least including tools.tsv & default.tsv & job configuration.tsv three types of configuration files__ |
| __data__ | Bundled data needed by the pipeline, such as genomes, indices, bcsubset and other software |
| __modules__ | Modular .smk scripts for the pipeline; the scripts folder under it contains tool scripts (py, R) | 
| logs | Storage for log files; it is recommended to use cmd &> logs/task.log to save logs here |
| inputs | Storage path for example input fq.gz files; your data can also be in other paths, __ignore__ |
| outputs | Storage path for example output folders; your output can also be in other paths, __ignore__ |

## Configs

Configuration file writing instructions

---

tools.tsv [Tool file configuration items]
| column | meaning | others |
| --- | --- | --- |
| tool | tool name | must be consistent with the default "tools.tsv" file, cannot be changed |
| path | tool access path | configure according to your environment; when setting up/migrating the environment, this usually needs to be reconfigured |

---

job.tsv & default.tsv [Task parameter configuration items]

_Note: The role of default.tsv is to fill in default values only when an item is missing in job.tsv; only the items listed below that can be configured with default parameters will take effect in default.tsv_

| column | meaning | can be configured as default parameter(default.tsv) |
| --- | --- | --- |
| Job | task name, also the name of the output files |  |
| DataDir | root directory of input files, for convenience of configuration, can be empty | √ |
| R1 | R1 path [barcode], if DataDir is configured, it becomes "{DataDir}/{R1}"	 |  |
| R2 | R2 path [read], if DataDir is configured, it becomes "{DataDir}/{R2}" |  |
| RU | RU path [UMI], if DataDir is configured, it becomes "{DataDir}/{RU}" |  |
| Prefix | __if R1,R2,RU are configured, this is invalid__, reads naming prefix, used to determine the input path "{DataDir}/{Prefix}_[R1/R2/RU].{Postfix}" for 2/3 kinds of inputs |  |
| Postfix | __if R1,R2,RU are configured, this is invalid__, reads naming suffix, used to determine the input path "{DataDir}/{Prefix}_[R1/R2/RU].{Postfix}" for 2/3 kinds of inputs | √ |
| Platform | data sequencing platform, supported parameters: 10x3-v1,10x3-v1.r3,10x3-v2,10x3-v3,10x5,10x5-v1,10x5-v2 | √ |
| Tissue | tissue, used for celltypist annotation, supports multiple tissues, e.g. [bone_marrow, pbmc] or [bone_marrow pbmc] separated by "," or " ", supported parameters are the tissue dictionary key names inside celltypist | √ |
| RefGenome | reference genome path (directory) | √ |
| RefGenomeIndex | reference genome index path (file) | √ |
| CellTypistModel | celltypist annotation model path (directory is enough) | √ |
| OutDir | output directory for this task, the output folder will be at path "{OutDir}/{Job}" | √ |
| Memory | memory parameter, passed to software that supports it | √ |
| Thread | core/thread parameter, passed to software that supports it | √ |
| CondaEnv | conda environment directory supporting sceasy's Seurat->H5ad conversion | √ |
| Description | description of this task, usually can be empty | √ |


## Outputs

Description of output directory items

| dir | description |
| --- | --- |
| 0.preFq | Due to the existence of v1 files, a unified fq preprocessing step is performed; if no special processing is needed, the pipeline will use ln -s to create soft links, __so be careful not to delete the source fq.gz, this output is not a full copy__, for v1 data an isalign.flag will be output to mark whether the seqID is aligned |
| 1.Map | STAR mapping output path |
| 2.Filter | Quality control filtering output path for the STAR output gene expression matrix, extracting high-quality cell barcodes, and also includes celltypist annotation results |
| 3.scTE | scTE output csv.gz |
| benchmark | Detailed time and CPU usage for each step, can be checked if resource requests need to be optimized |
| log | logs, including out and err, recording stdout and stderr |


