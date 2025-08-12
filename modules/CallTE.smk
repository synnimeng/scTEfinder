#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# # # # # # # # # # # # 
"""
╭───────────────────────────────────────╮ 
│ CallTE.smk   2025/03/03/-12:35 
╰───────────────────────────────────────╯ 
│ Description:
    输入bam, 进行 call TE, 生成TE数据
""" # [By: HuYw]

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
        expectCellN = 20000,    # 期待获取细胞量
    message:
        """
        CallTE by scTE for Job:{wildcards.job} {params.platform}.
        """
    shell:
        """
        cd {params.outdir}  # 必须 cd, 否则会输出至当前文件夹
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

