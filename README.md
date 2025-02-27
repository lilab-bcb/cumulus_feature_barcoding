# Cumulus Tool on Feature Barcoding

[![](https://img.shields.io/github/v/release/lilab-bcb/cumulus_feature_barcoding.svg)](https://github.com/lilab-bcb/cumulus_feature_barcoding/releases)

A fast C++ tool to extract feature-count matrix from sequence reads in FASTQ files. We uses isal-l for decompressing and Heng Li's kseq library for read parsing. It is used by Cumulus for feature-count matrix generation of cell hashing, nucleus hashing, CITE-Seq and Perturb-seq protocols, using either 10x Genomics V2, V3 or V4 chemistry.

## Installation

The installation has been tested on Debian and Ubuntu Linux.

1. Install dependency packages:

```
sudo apt install build-essential git libisal2 libisal-dev libdeflate0 libdeflate-dev
```

**Important:** Make sure to install `libisal2` and `libisal-dev` version 2.30.0 or later.

2. Check out this repository via Git:

```
git clone https://github.com/lilab-bcb/cumulus_feature_barcoding.git
```

3. Enter the directory and compile:

```
cd cumulus_feature_barcoding
make all
```

4. Now you'll have an executable named ``generate_count_matrix_ADTs`` inside your folder. Type

```
./generate_count_matrix_ADTs
```

to see its usage.

5. The 10x barcode inclusion list files are needed. Below is a list of them which you can find in Cell Ranger v9.0.1 source code:

* `cellranger-9.0.1/lib/python/cellranger/barcodes/`:
  * `3M-3pgex-may-2023_TRU.txt.gz`
  * `3M-5pgex-jan-2023.txt.gz`
  * `3M-february-2018_TRU.txt.gz`
  * `737K-arc-v1.txt.gz`
  * `737K-august-2016.txt`
* `cellranger-9.0.1/lib/python/cellranger/barcodes/translation/`:
  * `3M-3pgex-may-2023_NXT.txt.gz`
  * `3M-february-2018_NXT.txt.gz`

Just put these 7 files in the same folder, and specify it as `cell_barcodes_dir` argument for `generate_count_matrix_ADTs` command.

### Compile on Mac OS

1. Homebrew installation

```
brew install isa-l
brew install libdeflate
```

2. Compile

```
cd cumulus_feature_barcoding
make all
```

## Usage
