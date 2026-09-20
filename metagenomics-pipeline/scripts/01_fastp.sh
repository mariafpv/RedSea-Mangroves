#!/bin/bash

module load fastp 

source config/config.sh

COMPARTMENT="$1"

INPUT_DIR="${DATA_DIR}/${COMPARTMENT}"
OUTPUT_DIR="PreProcess/${COMPARTMENT}"

mkdir -p "${OUTPUT_DIR}"

for R1 in "${INPUT_DIR}"/*_R1_001.fastq.gz; do

    SAMPLE=$(basename "${R1}" _R1_001.fastq.gz)
    R2="${INPUT_DIR}/${SAMPLE}_R2_001.fastq.gz"

    if [[ ! -f "${R2}" ]]; then
        echo "Warning: Missing R2 for ${SAMPLE}"
        continue
    fi

    fastp \
        -i "${R1}" \
        -I "${R2}" \
        -o "${OUTPUT_DIR}/${SAMPLE}.filt_R1.fastq.gz" \
        -O "${OUTPUT_DIR}/${SAMPLE}.filt_R2.fastq.gz" \
        -e 30 \
        -w "${THREADS}" \
        -j "${OUTPUT_DIR}/${SAMPLE}.fastp.json" \
        -h "${OUTPUT_DIR}/${SAMPLE}.fastp.html"

done
