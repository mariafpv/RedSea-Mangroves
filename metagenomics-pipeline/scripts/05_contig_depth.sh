#!/bin/bash

source config/config.sh

COMPARTMENT="$1"

BAMDIR="${MAPPING_DIR}/${COMPARTMENT}_bam"
OUTDIR="${BINNING_DIR}/${COMPARTMENT}"

mkdir -p "${OUTDIR}"

jgi_summarize_bam_contig_depths \
    --outputDepth "${OUTDIR}/${COMPARTMENT}_bin_depth.txt" \
    --noIntraDepthVariance \
    "${BAMDIR}"/*.bam
