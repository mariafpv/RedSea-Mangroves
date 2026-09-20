#!/bin/bash

module load metawrap/1.3
module load pigz

source config/config.sh

COMPARTMENT="$1"

BINNING="${BINNING_DIR}/${COMPARTMENT}"
OUTDIR="${BINNING}/metawrap_refinement"
STATDIR="${BINNING}/stats"

mkdir -p "${OUTDIR}" "${STATDIR}"

metawrap bin_refinement \
    -t "${THREADS}" \
    -o "${OUTDIR}" \
    -A "${BINNING}/concoct/concoct_bins" \
    -B "${BINNING}/maxbin_bins" \
    -C "${BINNING}/metabat_bins" \
    -c 50 \
    -x 10

# Save MAG quality statistics
cp \
    "${OUTDIR}/metawrap_50_10_bins.stats" \
    "${STATDIR}/${COMPARTMENT}_MAG_stats.tsv"

# Rename MAGs
for BIN in "${OUTDIR}/metawrap_50_10_bins"/*.fa; do

    [[ -e "${BIN}" ]] || continue

    NAME=$(basename "${BIN}")

    mv \
        "${BIN}" \
        "${OUTDIR}/metawrap_50_10_bins/${COMPARTMENT}_${NAME}"

done

# Compress MAGs
pigz \
    -p "${THREADS}" \
    "${OUTDIR}/metawrap_50_10_bins/"*.fa
