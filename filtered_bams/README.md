# Mapping reads and filtering bams
## List of Scripts
####
#### Step 1: Map reads using Bowtie2 (paired-end for new data, single reads for reference bee sequences downloaded from ncbi)
-map_reads.sh
-map_reads_single.sh

#### Step 2: Sort reads using samtools
-sort_bams.sh

#### Step 3: De-duplicate reads using PICARD MarkDuplicates and calculate BAQ
-dedup_baq_bams.sh
-dedup_baq_bams_mac.sh # same except path to picard

#### Download reference bee data
-download_fastq_ena.sh

## Record of scripts run
- commands.txt # includes preliminary analyses not in the paper
