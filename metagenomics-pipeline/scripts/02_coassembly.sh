#!/bin/bash

source config/config.sh

COMPARTMENT="$1"

DATA="${DATA_DIR}/${COMPARTMENT}"
OUT="${COASSEMBLY_DIR}/${COMPARTMENT}"
TMP="${COASSEMBLY_DIR}/tmp/${COMPARTMENT}"

mkdir -p "${OUT}" "${TMP}"

R1_LIST=()
R2_LIST=()

for R1 in "${DATA}"/*.filt_R1.fastq.gz; do

    SAMPLE=$(basename "${R1}" .filt_R1.fastq.gz)
    R2="${DATA}/${SAMPLE}.filt_R2.fastq.gz"

    if [[ -f "${R2}" ]]; then
        R1_LIST+=("${R1}")
        R2_LIST+=("${R2}")
    else
        echo "Warning: Missing R2 for ${SAMPLE}"
    fi

done

if [[ ${#R1_LIST[@]} -eq 0 ]]; then
    echo "Error: No paired-end samples found."
    exit 1
fi

R1=$(IFS=,; echo "${R1_LIST[*]}")
R2=$(IFS=,; echo "${R2_LIST[*]}")

megahit \
    -1 "${R1}" \
    -2 "${R2}" \
    -o "${OUT}" \
    --tmp-dir "${TMP}" \
    --num-cpu-threads "${THREADS}" \
    --min-contig-len 1000
