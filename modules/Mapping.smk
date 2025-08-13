#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# # # # # # # # # # # # 
"""
╭───────────────────────────────────────╮ 
│ Mapping.smk   2025/03/01/-17:18 
╰───────────────────────────────────────╯ 
│ Description:
    Map paired FASTQ files from various platforms to the reference genome
""" # [By: Scripts-synni Meng, SMK-yiwen Hu]

# region |- Import -|
from types import SimpleNamespace
import configParser
import os
# endregion

# T = configParser.EnvConfig(config['tool'])
# Q = configParser.JobQueue(config['job'])

rule preprocessFastq:
    input: 
        R1 = lambda wildcards: Q.getjob(wildcards.job).R1FILE,
        R2 = lambda wildcards: Q.getjob(wildcards.job).R2FILE,
    output: 
        R1 = "{outdir}/{job}/0.preFq/{job}.R1.fastq.gz",
        R2 = "{outdir}/{job}/0.preFq/{job}.R2.fastq.gz",
    log:
        out = "{outdir}/{job}/log/{job}.0.preFq.out", 
        err = "{outdir}/{job}/log/{job}.0.preFq.err", 
    benchmark:
        "{outdir}/{job}/benchmark/{job}.0.preFq.benchmark",
    params:
        platform = lambda wildcards: Q.getjob(wildcards.job).Platform,
        thread = lambda wildcards: Q.getjob(wildcards.job).Thread,
        umiR = lambda wildcards: Q.getjob(wildcards.job).RUFILE,
    message:
        """
        Preprocess fastq for Job:{wildcards.job} {params.platform}.
        """
    shell:
        """
        if [ {params.platform} == "10x3-v1" ]; then
            bash modules/scripts/convert10xv1.sh \
                -B {input.R1} \
                -R {input.R2} \
                -X {output.R1} \
                -O {output.R2} \
                -t {params.thread} \
                1> {log.out} 2> {log.err}
            # {T.python} modules/scripts/10xv1Convert.py \
            #     -r1 {input.R1} \
            #     -r2 {input.R2} \
            #     -o1 {output.R1} \
            #     -o2 {output.R2} \
            #     -n {params.thread} \
            #     1> {log.out} 2> {log.err}
        elif [ {params.platform} == "10x3-v1.r3" ]; then
            bash modules/scripts/convert10xv1.3r.sh \
                -B {input.R1} \
                -R {input.R2} \
                -U {params.umiR} \
                -X {output.R1} \
                -O {output.R2} \
                -t {params.thread} \
                1> {log.out} 2> {log.err}
        else
            ln -sf {input.R1} {output.R1}
            ln -sf {input.R2} {output.R2}
        fi
        """



def STAR_ParamByPlatform(platform):
    platform = platform.split('.')[0]
    # whitelist
    WLDir = "data/whitelist"
    ParamsWL = {
        "10x3-v1": "737K-april-2014_rc.txt",
        "10x3-v2": "737K-august-2016.txt",
        "10x3-v3": "3M-february-2018.txt",
    }
    ParamsWL["10x5-v1"] = ParamsWL["10x5-v2"] = ParamsWL["10x5"] = ParamsWL["10x3-v2"]
    assert platform in ParamsWL.keys(), \
    f"Error: Job platform only support in {tuple(params.keys())}, but get {platform}"
    # barcode length param
    ParamsCB = {
        "10x3-v1": "--soloCBlen 14 --soloUMIstart 15 --soloUMIlen 10",
        "10x3-v2": "--soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10",
        "10x3-v3": "--soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12",
        #"10x5": "--soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10",
    }
    ParamsCB["10x5-v1"] = ParamsCB["10x5-v2"] = ParamsCB["10x5"] = ParamsCB["10x3-v2"]
    # other param
    ParamOther = ""
    if platform.startswith("10x3"):
        ParamOther+="--clipAdapterType CellRanger4"
    elif platform.startswith("10x5"):
        ParamOther+="--soloBarcodeMate 1 --clip5pNbases 39 0"
    # all params
    return f"{ParamsCB[platform]} {ParamOther} --soloCBwhitelist {WLDir}/{ParamsWL[platform]}"
    



rule MapBySTAR:
    input: 
        R1 = "{outdir}/{job}/0.preFq/{job}.R1.fastq.gz",
        R2 = "{outdir}/{job}/0.preFq/{job}.R2.fastq.gz",
    output: 
        tmp = temp(directory("{outdir}/{job}/1.Map/{job}._STARtmp")),
        Bam = "{outdir}/{job}/1.Map/{job}.Aligned.sortedByCoord.out.bam",
        Matrix = directory("{outdir}/{job}/1.Map/{job}.Solo.out/GeneFull_Ex50pAS/raw/"),
        FilterMatrix = directory("{outdir}/{job}/1.Map/{job}.Solo.out/GeneFull_Ex50pAS/filtered/"),
        # barcodes = "{outdir}/{job}/1.Map/{job}.Solo.out/GeneFull_Ex50pAS/filtered/barcodes.tsv",
    log:
        out = "{outdir}/{job}/log/{job}.1.Map.out", 
        err = "{outdir}/{job}/log/{job}.1.Map.err", 
    benchmark:
        "{outdir}/{job}/benchmark/{job}.1.Map.benchmark",
    params:
        starPrefix = "{outdir}/{job}/1.Map/{job}.",
        refGenome = lambda wildcards: Q.getjob(wildcards.job).RefGenome,    # genome
        platform = lambda wildcards: Q.getjob(wildcards.job).Platform,
        thread = lambda wildcards: Q.getjob(wildcards.job).Thread,
        memory = lambda wildcards: Q.getjob(wildcards.job).Memory,
        starParam = lambda wildcards: STAR_ParamByPlatform(Q.getjob(wildcards.job).Platform),
    message:
        """
        Mapping by STAR for Job:{wildcards.job} {params.platform}.
        """
    shell:
        """
        ulimit -n 10000
#        {T.STAR} --runThreadN {params.thread} 
        {T.STAR} --runThreadN 12 \
        --genomeDir {params.refGenome} \
        --readFilesCommand zcat --readFilesIn {input.R2} {input.R1} \
        --outFileNamePrefix {params.starPrefix} \
        --runDirPerm All_RWX \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterScoreMin 30 \
        --soloFeatures GeneFull_Ex50pAS Velocyto \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
        --soloCellFilter EmptyDrops_CR \
        --soloMultiMappers EM \
        --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --outMultimapperOrder Random --runRNGseed 777 --outSAMmultNmax 1 \
        --outSAMattributes NH HI AS nM CR CY UR UY CB UB \
        --soloType CB_UMI_Simple \
        --soloBarcodeReadLength 0 \
        --soloCBstart 1 {params.starParam} 1> {log.out} 2> {log.err}
        """


