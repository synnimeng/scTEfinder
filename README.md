[![DOI](https://zenodo.org/badge/1010980213.svg)](https://doi.org/10.5281/zenodo.17396260)


# scTEfinder 

scTEfinder processes data from FASTQ files through three core modules: 
* Read mapping,
* Quality control (QC),
* TE quantification.

Specifically, it aligns reads to the reference genome, keeps multi-mapped reads to capture TEs, filters BAM files based on cell QC, and counts TEs at the subfamily level. 
scTEfinder outputs a combined gene-TE count matrix that can be directly used in standard downstream analyses in Seurat (R) or Scanpy (Python).

![Pipeline](https://github.com/synnimeng/scTEfinder/blob/main/imgs/pipeline.png "scTEfinder Pipeline") 

Here, we present the snakemake implementation of scTEfinder for TE quantification from a set of fq.gz reads.
<br/>

### Environment

Use the following commands to build the scTEfinder environment with Anaconda.
For full details, see: [scTE-Conda.sh](https://github.com/synnimeng/scTEfinder/blob/main/scTE-Conda.sh)

```shell
# create conda env
conda create -n scTEfinder 'python=3.9' \
 bioconda::scrublet \
 bioconda::seqkit bioconda::star \
 bioconda::bedtools bioconda::samtools \
 conda-forge::r-seurat conda-forge::r-devtools \
 conda-forge::r-haven conda-forge::r-tidyverse conda-forge::r-argparse -y

# python libs
conda activate scTEfinder
TE_ENV=$(conda info --envs | grep '*' | awk '{print $3}')
export LD_LIBRARY_PATH="${TE_ENV}/lib:$LD_LIBRARY_PATH"
pip install scte
pip install scanpy
pip install scrublet
pip install cython
pip install snakemake==6 # / pip install snakemake==8
pip install --upgrade pulp==2.7.0
pip install celltypist

# Usually, numpy > 2 conflicts with scanpy—upgrade/downgrade as needed
pip install --upgrade 'numpy<2'

```

Additional R packages required:

```R
devtools::install_github("cellgeni/sceasy")
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
```
<br/>

### Running the pipeline

Full execution instructions are in: [runSMK.sh](https://github.com/synnimeng/scTEfinder/blob/main/runSMK.sh)

```shell
conda activate scTEfinder
# cd /path/to/scTEfinder

# run demo pipeline
snakemake --config job=configs/demo.tsv -j 4
```
<br/>

### Input, output, and configuration file overview

Please refer to the documentation: [MANUAL.md](https://github.com/synnimeng/scTEfinder/blob/main/MANUAL.md)

<br/>

### Additional notes

Large files in inputs and data are stored on Google Drive; please download and extract them locally with:

```shell
tar -zxf *.tar.gz
```















