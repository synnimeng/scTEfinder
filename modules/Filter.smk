#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# # # # # # # # # # # # 
"""
╭───────────────────────────────────────╮ 
│ Filter.smk   2025/03/01/-17:17 
╰───────────────────────────────────────╯ 
│ Description:
    输入mapping后路径, 作 bam & gene expr 质控过滤
""" # [By: HuYw]

# region |- Import -|
from types import SimpleNamespace
import configParser
import os
# endregion

# T = configParser.EnvConfig(config['tool'])
# Q = configParser.JobQueue(config['job'])
ruleorder: SeuratFilter > ScanpyFilter  # 优先使用Seurat
rule SeuratFilter:
    input: 
        STARsolo = "{outdir}/{job}/1.Map/{job}.Solo.out/GeneFull_Ex50pAS/filtered",
    output: 
        SeuratObject = "{outdir}/{job}/2.Filter/{job}.rds",
        Barcodes = "{outdir}/{job}/2.Filter/{job}.cells.csv",
        AnnData = "{outdir}/{job}/2.Filter/{job}.h5ad",
    log:
        out = "{outdir}/{job}/log/{job}.2.Filter.1.Seurat.out", 
        err = "{outdir}/{job}/log/{job}.2.Filter.1.Seurat.err", 
    benchmark:
        "{outdir}/{job}/benchmark/{job}.2.Filter.1.Seurat.benchmark",
    params:
        job = "{job}",
        outdir = "{outdir}/{job}/2.Filter/",
        platform = lambda wildcards: Q.getjob(wildcards.job).Platform,
        thread = lambda wildcards: Q.getjob(wildcards.job).Thread,
        condaEnv = lambda wildcards: Q.getjob(wildcards.job).CondaEnv,
    message:
        """
        Filter Cell by Seurat[R] for Job:{wildcards.job} {params.platform}.
        """
    shell:
        """
        {T.Rscript} modules/scripts/SeuratFilter.R -i {input.STARsolo} -s {params.job} \
        -c {params.thread} -o {params.outdir} -e {params.condaEnv} \
        1> {log.out} 2> {log.err}
        # -p {params.platform} 
        """



rule ScanpyFilter:
    input: 
        STARsolo = "{outdir}/{job}/1.Map/{job}.Solo.out/GeneFull_Ex50pAS/filtered",
    output: 
        SeuratObject = "{outdir}/{job}/2.Filter/{job}.h5ad",
        Barcodes = "{outdir}/{job}/2.Filter/{job}.cells.csv",
    log:
        out = "{outdir}/{job}/log/{job}.2.Filter.1.Scanpy.out", 
        err = "{outdir}/{job}/log/{job}.2.Filter.1.Scanpy.err", 
    benchmark:
        "{outdir}/{job}/benchmark/{job}.2.Filter.1.Scanpy.benchmark",
    params:
        job = "{job}",
        outdir = "{outdir}/{job}/2.Filter/",
        platform = lambda wildcards: Q.getjob(wildcards.job).Platform,
        thread = lambda wildcards: Q.getjob(wildcards.job).Thread,
    message:
        """
        Filter Cell by Scanpy[py] for Job:{wildcards.job} {params.platform}.
        """
    shell:
        """
        {T.python} modules/scripts/ScanpyFilter.py -i {input.STARsolo} -j {params.job} \
        -c {params.thread} -o {params.outdir} \
        1> {log.out} 2> {log.err}
        # -p {params.platform} 
        """


rule SubsetBam:
    input: 
        Bam = "{outdir}/{job}/1.Map/{job}.Aligned.sortedByCoord.out.bam",
        Barcodes = "{outdir}/{job}/2.Filter/{job}.cells.csv",
    output: 
        Bam = "{outdir}/{job}/2.Filter/{job}.sub.bam",
    log:
        out = "{outdir}/{job}/log/{job}.2.Filter.2.Subset.out", 
        err = "{outdir}/{job}/log/{job}.2.Filter.2.Subset.err", 
    benchmark:
        "{outdir}/{job}/benchmark/{job}.2.Filter.2.Subset.benchmark",
    params:
        platform = lambda wildcards: Q.getjob(wildcards.job).Platform,
        thread = lambda wildcards: Q.getjob(wildcards.job).Thread,
    message:
        """
        Subset Bam in filter cells for Job:{wildcards.job} {params.platform}.
        """
    shell:
        """
        {T.bcsubset} -w {input.Barcodes} -o {output.Bam} {input.Bam} 1> {log.out} 2> {log.err}
        """
