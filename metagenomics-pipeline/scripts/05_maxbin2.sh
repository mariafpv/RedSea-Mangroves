#!/bin/bash

source config/config.sh

COMPARTMENT="$1"

CONTIG="${MAPPING_DIR}/${COMPARTMENT}/${COMPARTMENT}_contigs.fa"
DEPTH="${BINNING_DIR}/${COMPARTMENT}/${COMPARTMENT}_bin_depth.txt"
OUTDIR="${BINNING_DIR}/${COMPARTMENT}/maxbin_bins"

mkdir -p "${OUTDIR}"

run_MaxBin.pl \
    -contig "${CONTIG}" \
    -abund "${DEPTH}" \
    -out "${OUTDIR}/${COMPARTMENT}_maxbin" \
    -thread "${THREADS}"

for BIN in "${OUTDIR}/${COMPARTMENT}_maxbin".*.fasta; do

    [[ -e "${BIN}" ]] || continue

    NAME=$(basename "${BIN}")

    mv \
        "${BIN}" \
        "${OUTDIR}/${COMPARTMENT}_${NAME}"

done
