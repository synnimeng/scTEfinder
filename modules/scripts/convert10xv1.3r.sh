#!/bin/sh
<<README
[convert10x_v1.sh] => 2025/03/01/-17:56	 # by synnimeng & yiwen hu
Intro:	convert 10x 3' v1 data.
Usage:	bash convert10x_v1.sh -B <barcode.fq.gz> -R <read.fq.gz> -U <UMI.fq.gz> -X <barcode+umi/r1> -O <read/r2>"
        # B: barcode, R: read, U: UMI

README


# Define help function
usage() {
    echo "Usage: $0 -r1 <R1.fq.gz> -r2 <R2.fq.gz> -i1 <I1.fq.gz> -o <output_dir>"
    echo "Options:"
    echo "  -B      Input Barcode reads file (fastq.gz)"
    echo "  -R      Input reads file (fastq.gz)"
    echo "  -U      Input UMI reads file (fastq.gz)"
    echo "  -X      Output R1 file"
    echo "  -O      Output R2 file"
    echo "  -t       [option] Max Cpu Threads to use, 4"
    exit 1
}

# Parse command-line arguments
while getopts ":B:R:U:X:O:t::" opt; do
    echo $opt $OPTARG
    case ${opt} in
        B ) inR1=$OPTARG ;;
        R ) inR2=$OPTARG ;;
        U ) inI1=$OPTARG ;;
        X ) oR1=$OPTARG ;;
        O ) oR2=$OPTARG ;;
        t  ) nThread=${OPTARG:-4} ;;
        \? ) usage ;;
    esac
done


# Check if required arguments are provided
if [[ -z "$inR1" || -z "$inR2" || -z "$inI1" || -z "$oR1" || -z "$oR2" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Create output directory
mkdir -p "$(dirname $oR1)"
mkdir -p "$(dirname $oR2)"
# Temporary folder
tmpDir=$(dirname $oR1)
tmpSID1="$tmpDir/SID1.txt"
tmpSID2="$tmpDir/SID2.txt"
tmpUMI="$tmpDir/UMI.fastq.gz"
tmpAlignFlag="$tmpDir/isalign.flag"

# Extract UMI and verify read sequence IDs
echo "10X 3' v1 prepare handling..."
seqkit seq -j $nThread "$inR1" -n -i > "$tmpSID1" &
seqkit seq -j $nThread "$inR2" -n -i > "$tmpSID2" &
zcat "$inI1" | seqkit subseq -j $nThread -r 1:10 -o "$tmpUMI"   # Take first 10 bp as UMI sequence

wait    # Wait for the two child SeqID counting processes to finish
# Compare whether IDs are identical
if diff "$tmpSID1" "$tmpSID2" > /dev/null; then
    echo "IDs are identical."
    echo "True" > "$tmpAlignFlag"   # Whether IDs are aligned
else
    echo "Warning: IDs are different!"
    echo "False" > "$tmpAlignFlag"  # 
    exit 1
fi

rm "$tmpSID1" "$tmpSID2"

#  Whether IDs are aligned
echo "Concatenating CB and UMI files..."
zcat "$inR1" | paste - - - - <(zcat "$tmpUMI" | paste - - - -) | \
    awk 'BEGIN{OFS="\n"} {print $1" "$2, $3 $8, $4, $5 $10}' | \
    gzip -c > "$oR1"

# Force overwrite and create symbolic link for R2 file
ln -sf "$inR2" "$oR2"

rm "$tmpUMI"

echo "Processing complete."
