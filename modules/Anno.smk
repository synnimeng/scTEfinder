#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# # # # # # # # # # # # 
"""
╭───────────────────────────────────────╮ 
│ Anno.smk   2025/03/04/-17:37 
╰───────────────────────────────────────╯ 
│ Description:
    对 h5ad AnnData 执行 细胞类型 注释
""" # [By: HuYw]

# region |- Import -|
from types import SimpleNamespace
import configParser
import os
# endregion

# T = configParser.EnvConfig(config['tool'])
# Q = configParser.JobQueue(config['job'])

rule AnnoByCellTypist:
    input: 
        AnnData = "{outdir}/{job}/2.Filter/{job}.h5ad",
    output: 
        Label = "{outdir}/{job}/2.Filter/{job}.pred_labels.csv",
    log:
        out = "{outdir}/{job}/log/{job}.2.Filter.3.Anno.out", 
        err = "{outdir}/{job}/log/{job}.2.Filter.3.Anno.err", 
    benchmark:
        "{outdir}/{job}/benchmark/{job}.2.Filter.3.Anno.benchmark",
    params:
        prefix = "{outdir}/{job}/2.Filter/{job}",
        model = lambda wildcards: Q.getjob(wildcards.job).CellTypistModel,
        tissue = lambda wildcards: Q.getjob(wildcards.job).Tissue,
        platform = lambda wildcards: Q.getjob(wildcards.job).Platform,
    message:
        """
        Anno Cell by CellTypist[py] for Job:{wildcards.job} {params.platform}.
        """
    shell:
        """
        {T.python} modules/scripts/runCellTypist.py -i {input.AnnData} -m {params.model} \
        -t {params.tissue} -o {params.prefix} -d \
        1> {log.out} 2> {log.err}
        """