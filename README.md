# scTEfinder 
__(snakemake pipeline)__

This pipeline is designed to call TEs from a set of fq.gz reads.

请在runSMK.sh 脚本查看运行步骤
而scTE-Conda.sh 提供了运行时环境的完整搭建步骤

## Description

本项目目录功能说明

| dir | description |
| --- | --- |
| __configs__ | 流程所需的配置文件, __至少包含 tools.tsv & default.tsv & 任务配置.tsv__ 3种配置文件 |
| __data__ | 流程所需的绑定数据, 如基因组、index、bcsubset等软件 |
| __modules__ | 流程的模块化.smk脚本, 其下scripts文件夹包含工具脚本(py、R) | 
| logs | 日志文件储存地, 推荐使用 cmd &> logs/task.log 将日志储存至此 |
| inputs | 示例输入fq.gz文件的储存路径, 您的数据也可以在其他路径, __忽略此项__ |
| outputs | 示例输出文件夹的储存路径，您的输出也可以在其他路径，__忽略此项__ |

## Configs

配置文件编写说明

---

tools.tsv [工具文件配置项]
| 列 | 释义 | 其他 |
| --- | --- | --- |
| tool | 工具名称 | 与默认文件提供一致, 不可更改 |
| path | 工具访问路径 | 根据个人环境进行配置, 当进行环境配置/迁移, 一般需要重新配置 |

---

job.tsv & default.tsv [任务参数配置项]

_注: default.tsv的作用是 仅在job.tsv缺失某项信息时进行默认值填充; 仅下方 可配置默认参数的项目, 才可在default.tsv生效_
| 列 | 释义 | 是否可配置默认参数 |
| --- | --- | --- |
| Job | 任务名称, 也是 输出文件的命名 |  |
| DataDir | 输入文件的根目录, 仅为了方便配置, 可为空 | √ |
| R1 | R1路径[barcode], 若配置DataDir, 为"{DataDir}/{R1}" |  |
| R2 | R2路径[read], 若配置DataDir, 为"{DataDir}/{RW}" |  |
| RU | RU路径[UMI], 若配置DataDir, 为"{DataDir}/{RU}" |  |
| Prefix | __若配置R1,R2,RU, 则失效__, reads的命名前缀, 将据此确定输入路径"{DataDir}/{Prefix}_[R1/R2/RU].{Postfix}" 2/3种输入 |  |
| Postfix | __若配置R1,R2,RU, 则失效__, reads的命名后缀, 将据此确定输入路径"{DataDir}/{Prefix}_[R1/R2/RU].{Postfix}" 2/3种输入 | √ |
| Platform | 数据测序平台, 支持的参数有: 10x3-v1,10x3-v1.r3,10x3-v2,10x3-v3,10x5,10x5-v1,10x5-v2 | √ |
| Tissue | 组织, 用于celltypist标注, 支持多tissue传入, 如:[bone_marrow, pbmc]或[bone_marrow pbmc]使用","或" "分隔均可, 支持的参数是celltypist内部的tissue字典键名称 | √ |
| RefGenome | 参考基因组路径(目录) | √ |
| RefGenomeIndex | 参考基因组索引路径(文件) | √ |
| CellTypistModel | celltypist标注模型路径(目录即可) | √ |
| OutDir | 该任务输出目录, 输出文件夹将在路径"{OutDir}/{Job}" | √ |
| Memory | 若有软件支持 内存参数, 传入该值 | √ |
| Thread | 若有软件支持 核心/线程参数, 传入该值 | √ |
| CondaEnv | 支持sceasy的Seurat->H5ad 转换的conda环境目录 | √ |
| Description | 对该任务的描述, 一般为空 | √ |


## Outputs

输出目录项目说明


| dir | description |
| --- | --- |
| 0.preFq | 因v1文件的存在, 统一进行fq预处理操作, 若不需要特别处理, 流程将使用ln -s创建软链接, __因此谨慎删除源路径的fq.gz, 该输出并不是完全的copy__, 其中针对v1数据将输出_isalign.flag_标注seqID是否对齐 |
| 1.Map | STAR mapping的输出路径 |
| 2.Filter | 对STAR输出基因表达矩阵的质控过滤输出路径, 提取高质量细胞barcode, 也包含celltypist的标注结果 |
| 3.scTE | scTE输出 csv.gz |
| benchmark | 每个步骤详细的耗时及CPU占用, 若需要可以查看优化资源申请 |
| log | 日志, 包含out和err两种, 记录标准输出和错误 |

