#!/bin/bash

conda activate Samtools 

source config/config.sh

COMPARTMENT="$1"
SAMPLE="$2"

INDEX="${MAPPING_DIR}/${COMPARTMENT}/${COMPARTMENT}_contigs"
INPUT_DIR="PreProcess/${COMPARTMENT}"
OUTPUT_DIR="${MAPPING_DIR}/${COMPARTMENT}_bam"

mkdir -p "${OUTPUT_DIR}"

R1="${INPUT_DIR}/${SAMPLE}.filt_R1.fastq.gz"
R2="${INPUT_DIR}/${SAMPLE}.filt_R2.fastq.gz"

bowtie2 \
    --threads "${THREADS}" \
    -x "${INDEX}" \
    -1 "${R1}" \
    -2 "${R2}" \
    | samtools view -@ "${THREADS}" -b \
    | samtools sort \
        -@ "${THREADS}" \
        -o "${OUTPUT_DIR}/${SAMPLE}_sort.bam"

samtools index "${OUTPUT_DIR}/${SAMPLE}_sort.bam"
