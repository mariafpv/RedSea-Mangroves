#!/bin/bash

conda activate Samtools

source config/config.sh

COMPARTMENT="$1"

FASTA="${MAPPING_DIR}/${COMPARTMENT}/${COMPARTMENT}_contigs.fa"
INDEX="${MAPPING_DIR}/${COMPARTMENT}/${COMPARTMENT}_contigs"

mkdir -p "${MAPPING_DIR}/${COMPARTMENT}"

bowtie2-build \
    --threads "${THREADS}" \
    "${FASTA}" \
    "${INDEX}"
