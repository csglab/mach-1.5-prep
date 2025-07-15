# mach-1.5-prep

This repository contains scripts for downloading and preparing reference annotations and genomes from the NCBI RefSeq database. The processed data is used to train models in the [`mach-1.5-savanna`](https://github.com/your-org/mach-1.5-savanna) repository.

## Overview

- `scripts/`: Contains all preprocessing and utility scripts, including:
  - Downloading reference genomes and annotations (`download_references.py`)
  - Preparing datasets with or without UTRs (`prepare_refseq*.R`)
  - Generating phylogenetic codes, computing stats, and sampling subsets
  - Splitting data and converting it to `.jsonl.zst` format (`csv_gz_to_jsonl_zst.sh`)

- `data/refseq/`: Directory for storing downloaded and processed RefSeq data, organized into:
  - `annotation/`, `genome/`, `metadata/`, `sequences/`, etc.

## Output Format

The final output is a split dataset in `.jsonl.zst` compressed format compatible with the Savanna training framework.

---

For training with Mach-1.5, see [`mach-1.5-savanna`](https://github.com/goodarzilab/mach-1.5-savanna).

