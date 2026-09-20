#!/bin/bash

module load metabat/2.15.0

source config/config.sh

COMPARTMENT="$1"

CONTIG="${MAPPING_DIR}/${COMPARTMENT}/${COMPARTMENT}_contigs.fa"
DEPTH="${BINNING_DIR}/${COMPARTMENT}/${COMPARTMENT}_bin_depth.txt"
OUTDIR="${BINNING_DIR}/${COMPARTMENT}/metabat_bins"

mkdir -p "${OUTDIR}"

metabat2 \
    -i "${CONTIG}" \
    -a "${DEPTH}" \
    -o "${OUTDIR}/${COMPARTMENT}_metabat" \
    -t "${THREADS}"
