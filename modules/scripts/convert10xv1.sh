#!/bin/sh
<<README
[convert10x_v1.sh] => 2025/03/01/-17:56	 # by huyw
Intro:	来自Kimi生成, 参考xini姐代码, 转换10x 3' v1 数据
Usage:	bash convert10x_v1.sh -B <barcode.fq.gz> -R <read.fq.gz> -X <barcode+umi/r1> -O <read/r2>"
        # B: barcode, R: [read + umi]

README


# 定义帮助函数
usage() {
    echo "Usage: $0 -r1 <R1.fq.gz> -r2 <R2.fq.gz> -i1 <I1.fq.gz> -o <output_dir>"
    echo "Options:"
    echo "  -B      Input Barcode reads file (fastq.gz)"
    echo "  -R      Input reads file (fastq.gz)"
    echo "  -X      Output R1 file"
    echo "  -O      Output R2 file"
    echo "  -t       [option] Max Cpu Threads to use, 4"
    exit 1
}

# 解析命令行参数
while getopts ":B:R:X:O:t::" opt; do
    echo $opt $OPTARG
    case ${opt} in
        B ) inR1=$OPTARG ;;
        R ) inR2=$OPTARG ;;
        X ) oR1=$OPTARG ;;
        O ) oR2=$OPTARG ;;
        t  ) nThread=${OPTARG:-4} ;;
        \? ) usage ;;
    esac
done

<<EOF
inR1=ACGCGGAA_L001_R1.fastq.gz
inR2=ACGCGGAA_L001_R2.fastq.gz
oR1=R1.fq.gz
oR2=R2.fq.gz
nThread=4
EOF

# 检查必要的参数是否已提供
if [[ -z "$inR1" || -z "$inR2" || -z "$oR1" || -z "$oR2" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# 创建输出目录
mkdir -p "$(dirname $oR1)"
mkdir -p "$(dirname $oR2)"
# 临时文件夹
tmpDir=$(dirname $oR1)

# 定义输出文件路径
tmpSID1="$tmpDir/SID1.txt"
tmpSID2="$tmpDir/SID2.txt"
tmpUMI="$tmpDir/UMI.fastq.gz"
tmpAlignFlag="$tmpDir/isalign.flag"

# 提取 UMI 和校验 reads 序列ID号
echo "10X 3' v1 prepare handling..."
# 拆分R2 至 Read 和 UMI
zcat $inR2 | paste - - - - - - - - | \
    awk 'BEGIN{OFS="\n"} {print $1" "$2, $3, $4, $5}' | \
    gzip -c > "$oR2" &

zcat $inR2 | paste - - - - - - - - | \
    awk 'BEGIN{OFS="\n"} {print $6" "$7, $8, $9, $10}' | \
    gzip -c > "$tmpUMI" &


wait # 等待oR2 & UMI分发完毕
seqkit seq -j $nThread "$inR1" -n -i > "$tmpSID1" &
seqkit seq -j $nThread "$oR2" -n -i > "$tmpSID2" &


# 合并 Cell Barcode 和 UMI
echo "Concatenating CB and UMI files..."
zcat "$inR1" | paste - - - - <(zcat "$tmpUMI" | paste - - - -) | \
    awk 'BEGIN{OFS="\n"} {print $1" "$2, $3 $8, $4, $5 $10}' | \
    gzip -c > "$oR1"

wait    # 等待统计SeqID 
# 比较 ID 是否一致
if diff "$tmpSID1" "$tmpSID2" > /dev/null; then
    echo "IDs are identical."
    echo "True" > "$tmpAlignFlag"   # ID  是否对齐
else
    echo "Warning: IDs are different!"
    echo "False" > "$tmpAlignFlag"  # ID  是否对齐
    exit 1
fi

rm "$tmpSID1" "$tmpSID2"
rm "$tmpUMI"

echo "Processing complete."