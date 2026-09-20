#!/bin/bash

# Project directories
DATA_DIR="Data_combine"
MAPPING_DIR="Mapping"
BINNING_DIR="Binning"
COASSEMBLY_DIR="Coassembly"

# Computational resources
THREADS=${SLURM_CPUS_PER_TASK:-16}
