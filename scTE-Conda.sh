
conda create -n scTEfinder 'python=3.9' \
 hcc::aspera-cli bioconda::scrublet \
 bioconda::seqkit bioconda::star \
 bioconda::bedtools bioconda::samtools \
 conda-forge::r-seurat bioconda::r-archr \
 conda-forge::r-haven conda-forge::r-tidyverse conda-forge::r-argparse -y


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

# R
devtools::install_github("cellgeni/sceasy")
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')


# bc-subset [INSTALL TUTORIAL]
https://github.com/kehrlab/bcsubset/blob/main/README.md
# Error 1:  unrecognized command line option ‘-std=c++14’
# Error 2:  /lib64/libstdc++.so.6: version `GLIBCXX_3.4.30' not found (required by ./bcsubset)
<< EOF
conda activate scTEfinder 
# Error1
conda install -c conda-forge gcc_linux-64 gxx_linux-64
conda install zlib
cd $TE_ENV/bin
ln -s x86_64-conda-linux-gnu-gcc gcc
ln -s x86_64-conda-linux-gnu-g++ g++

# Error2
conda install -c anaconda libstdcxx-ng
# regist this to your ~/.bashrc, or as a hook fuction to call when you use 'conda activate te'
export LD_LIBRARY_PATH=${CONDA_HOME}/envs/te/lib:$LD_LIBRARY_PATH
EOF

# scTE [INSTALL TUTORIAL]
https://github.com/JiekaiLab/scTE
