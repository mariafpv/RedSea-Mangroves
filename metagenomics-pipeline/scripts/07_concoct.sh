#!/bin/bash

module load concoct
module load samtools

source config/config.sh

COMPARTMENT="$1"

CONTIG="${MAPPING_DIR}/${COMPARTMENT}/${COMPARTMENT}_contigs.fa"
BAMDIR="${MAPPING_DIR}/${COMPARTMENT}_bam"
OUTDIR="${BINNING_DIR}/${COMPARTMENT}/concoct"

mkdir -p "${OUTDIR}/concoct_bins"

# Fragment contigs
cut_up_fasta.py \
    "${CONTIG}" \
    -c 10000 \
    --merge_last \
    -b "${OUTDIR}/assembly_10K.bed" \
    -o 0 \
    > "${OUTDIR}/assembly_10K.fa"

# Calculate coverage
concoct_coverage_table.py \
    "${OUTDIR}/assembly_10K.bed" \
    "${BAMDIR}"/*.bam \
    > "${OUTDIR}/concoct_depth.txt"

# Run CONCOCT
concoct \
    -l 1500 \
    -t "${THREADS}" \
    --coverage_file "${OUTDIR}/concoct_depth.txt" \
    --composition_file "${OUTDIR}/assembly_10K.fa" \
    -b "${OUTDIR}/"

# Merge fragment clusters
merge_cutup_clustering.py \
    "${OUTDIR}/clustering_gt1500.csv" \
    > "${OUTDIR}/clustering_gt1500_merged.csv"

# Extract bins
python scripts/metawrap_split_concoct_bins.py \
    "${OUTDIR}/clustering_gt1500_merged.csv" \
    "${CONTIG}" \
    "${OUTDIR}/concoct_bins"
