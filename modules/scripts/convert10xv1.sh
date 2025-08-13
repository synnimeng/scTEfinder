#!/bin/sh
<<README
[convert10x_v1.sh] => 2025/03/01/-17:56	 # by synnimeng & yiwen hu
Intro:	convert 10x 3' v1 data.
Usage:	bash convert10x_v1.sh -B <barcode.fq.gz> -R <read.fq.gz> -X <barcode+umi/r1> -O <read/r2>"
        # B: barcode, R: [read + umi]

README


# Define help function
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

# Parse command-line arguments
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

# Check if required arguments are provided
if [[ -z "$inR1" || -z "$inR2" || -z "$oR1" || -z "$oR2" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Create output directory
mkdir -p "$(dirname $oR1)"
mkdir -p "$(dirname $oR2)"
# Temporary 
tmpDir=$(dirname $oR1)
tmpSID1="$tmpDir/SID1.txt"
tmpSID2="$tmpDir/SID2.txt"
tmpUMI="$tmpDir/UMI.fastq.gz"
tmpAlignFlag="$tmpDir/isalign.flag"

# Extract UMI and verify read sequence IDs
echo "10X 3' v1 prepare handling..."
# Split R2 into Read and UMI
zcat $inR2 | paste - - - - - - - - | \
    awk 'BEGIN{OFS="\n"} {print $1" "$2, $3, $4, $5}' | \
    gzip -c > "$oR2" &

zcat $inR2 | paste - - - - - - - - | \
    awk 'BEGIN{OFS="\n"} {print $6" "$7, $8, $9, $10}' | \
    gzip -c > "$tmpUMI" &


wait # 
seqkit seq -j $nThread "$inR1" -n -i > "$tmpSID1" &
seqkit seq -j $nThread "$oR2" -n -i > "$tmpSID2" &


# Concatenate Cell Barcode and UMI
echo "Concatenating CB and UMI files..."
zcat "$inR1" | paste - - - - <(zcat "$tmpUMI" | paste - - - -) | \
    awk 'BEGIN{OFS="\n"} {print $1" "$2, $3 $8, $4, $5 $10}' | \
    gzip -c > "$oR1"

wait 
# Compare whether IDs are identical
if diff "$tmpSID1" "$tmpSID2" > /dev/null; then
    echo "IDs are identical."
    echo "True" > "$tmpAlignFlag"   # Whether IDs are aligned
else
    echo "Warning: IDs are different!"
    echo "False" > "$tmpAlignFlag"  # Whether IDs are aligned
    exit 1
fi

rm "$tmpSID1" "$tmpSID2"
rm "$tmpUMI"

echo "Processing complete."
