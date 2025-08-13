#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# # # # # # # # # # # # 
"""
╭───────────────────────────────────────╮ 
│ Snakefile   2025/03/03/-01:14 
╰───────────────────────────────────────╯ 
│ Description:
    Whole pipeline from FASTQ → TE analysis, supports modular execution
""" # [By: ]

# region |- Import -|
from types import SimpleNamespace
import configParser
import os
# endregion

# region |- Global -|
def configset(config:dict, k:str, v):
    if k not in config:
        config[k] = v


assert 'job' in config.keys(), \
"Error: No config 'job' input, please use:\n'snakemake --config job=configs/job.tsv -j 2' to start a job!"
configset(config, 'tool', 'configs/tools.tsv')  # Tool configuration
configset(config, 'default', 'configs/default.tsv')  # Default-value configuration
# endregion


T = configParser.Tool(config['tool'])
Q = configParser.JobQueue(config['job'], config['default'])

# region |- Include .smk -|
include: "modules/Mapping.smk"
include: "modules/Filter.smk"
include: "modules/CallTE.smk"
include: "modules/Anno.smk"
# endregion


rule all:
    input:
        # # STAR map
        # expand("{outdir}/{job}/1.Map/{job}.Aligned.sortedByCoord.out.bam", zip, outdir=Q['OutDir'], job=Q['Job']),
        # # Filter Bam
        # expand("{outdir}/{job}/2.Filter/{job}.sub.bam", zip, outdir=Q['OutDir'], job=Q['Job']),
        # Anno Cell
        expand("{outdir}/{job}/2.Filter/{job}.pred_labels.csv", zip, outdir=Q['OutDir'], job=Q['Job']),
        # Call TE [scTE]
        expand("{outdir}/{job}/3.scTE/{job}.csv.gz", zip, outdir=Q['OutDir'], job=Q['Job']),
    message:
        "[Main] Whole Pipeline start..."


rule RunMapOnly:
    input:
        expand("{outdir}/{job}/1.Map/{job}.Aligned.sortedByCoord.out.bam", zip, outdir=Q['OutDir'], job=Q['Job']),
    message:
        "[Main] Run Map Only Pipeline start..."


rule RunFilterOnly:
    input:
        expand("{outdir}/{job}/2.Filter/{job}.sub.bam", zip, outdir=Q['OutDir'], job=Q['Job']),
    message:
        "[Main] Run Filter Only Pipeline start..."


rule RunAnnoOnly:
    input:
        expand("{outdir}/{job}/2.Filter/{job}.pred_labels.csv", zip, outdir=Q['OutDir'], job=Q['Job']),
    message:
        "[Main] Run Anno Only Pipeline start..."
