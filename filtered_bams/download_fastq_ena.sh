#!/bin/bash

# A script to download fastq files from European Nucleotide Archive
# right now only works with single-end reads (1 fastq)

# to run:
#filtered_bams$ ./download_fastq_ena.sh SRR5580833 ../data/Cridland_2018/fastq_files

# general bash script settings to make sure if any errors in the pipeline fail
# the it's a 'fail' and it passes all errors to exit and allows no unset variables
set -o pipefail
set -o errexit
set -o nounset

SRR="$1"
DIR="$2"
ftp_SRR=$(curl --silent "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$SRR&result=read_run&fields=fastq_ftp&download=txt" | tail -n 1)
echo "ftp file web address: $ftp_SRR"

aspera_SRR=$(echo "$ftp_SRR" | awk '{ gsub("ftp.sra.ebi.ac.uk", "era-fasp@fasp.sra.ebi.ac.uk:"); print }')
echo "aspera file web address: $aspera_SRR"


# download to output directory
mkdir -p "$DIR"
echo "downloading fastq to $DIR"
ascp -QT -l 300m -P33001 -i ~/Applications/Aspera\ Connect.app/Contents/Resources/asperaweb_id_dsa.openssh "$aspera_SRR" "$DIR"/.


# also download md5sum to check file integrity
md5=$(curl --silent "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$SRR&result=read_run&fields=fastq_md5,fastq_ftp&download=txt" | tail -n 1)
echo "md5sum: $md5 SRR: $SRR"
echo $md5 > "$DIR/$SRR.true.md5"
