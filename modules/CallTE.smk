#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# # # # # # # # # # # # 
"""
╭───────────────────────────────────────╮ 
│ CallTE.smk   2025/03/03/-12:35 
╰───────────────────────────────────────╯ 
│ Description:
    Input BAM file, call transposable elements (TEs) and generate TE data
""" # [By: Scripts-synni Meng, SMK-yiwen Hu]

# region |- Import -|
from types import SimpleNamespace
import configParser
import os
# endregion

# T = configParser.EnvConfig(config['tool'])
# Q = configParser.JobQueue(config['job'])


rule CallscTE:
    input: 
        Bam = "{outdir}/{job}/2.Filter/{job}.sub.bam",
    output: 
        TE = "{outdir}/{job}/3.scTE/{job}.csv.gz",
    log:
        out = "{outdir}/{job}/log/{job}.3.scTE.out", 
        err = "{outdir}/{job}/log/{job}.3.scTE.err", 
    benchmark:
        "{outdir}/{job}/benchmark/{job}.3.scTE.benchmark",
    params:
        outdir = "{outdir}/{job}/3.scTE/",
        TEprefix = "{job}",
        refGenomeIndex = lambda wildcards: Q.getjob(wildcards.job).RefGenomeIndex,
        platform = lambda wildcards: Q.getjob(wildcards.job).Platform,
        thread = lambda wildcards: Q.getjob(wildcards.job).Thread,
        expectCellN = 20000,    # Expected number of cells to recover
    message:
        """
        CallTE by scTE for Job:{wildcards.job} {params.platform}.
        """
    shell:
        """
        cd {params.outdir}  # do cd; otherwise output will be placed in the current directory
        {T.scTE} -p {params.thread} -x {params.refGenomeIndex} \
        -i {input.Bam} \
        -o {params.TEprefix} \
        -CB CB -UMI UB \
        --expect-cells {params.expectCellN} \
        --min_counts 0 \
        --min_genes 0 1> {log.out} 2> {log.err}
        # gzip
        gzip {params.TEprefix}.csv
        """

