
# app.conda
# conda activate te

cd /path/to/scTEfinder

TE_ENV='/path/to/envs/scte'
export LD_LIBRARY_PATH="${TE_ENV}/lib:$LD_LIBRARY_PATH"
export PATH="${TE_ENV}/bin:$PATH"

SMK='/path/to/snakemake'

# 查看流程图, 不实际运行
#$SMK --config job=configs/test.tsv -j 3 --dryrun

# 实际运行  # -j {3} 代表一次性只允许并行跑 {3} 个任务, 根据 资源申请 及 最大单任务耗费 而灵活变动
$SMK --config job=configs/test2.tsv -j 2 --unlock
$SMK --config job=configs/test2.tsv -j 2

# 查看所有支持的pipeline
# $SMK --config job=configs/demo.tsv --list-target-rules

#<<EOF
# 仅运行 Anno, 并刷新存在, 导致后续输出下次要 更新重跑
#$SMK --config job=configs/10x3v1.r3.tsv -j 3 --dryrun --allowed-rules RunAnnoOnly

# snakemake 只考虑文件 时间戳, 而不考虑代码更新进行rerun
#$SMK --config job=configs/10x3v1.r3.tsv -j 3 --dryrun --rerun-triggers mtime

# 修改输出 时间戳
#stat file
#touch -t 202503010113.49 file
#EOF
