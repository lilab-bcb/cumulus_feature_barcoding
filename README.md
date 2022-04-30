# Cumulus Tool on Feature Barcoding

[![](https://img.shields.io/github/v/release/lilab-bcb/cumulus_feature_barcoding.svg)](https://github.com/lilab-bcb/cumulus_feature_barcoding/releases)

A fast C++ tool to extract feature-count matrix from sequence reads in FASTQ files, which uses Heng Li's kseq library for read parsing. It is used by Cumulus for feature-count matrix generation of cell hashing, nucleus hashing, CITE-Seq and Perturb-seq protocols, using either 10x Genomics V2 or V3 chemistry.

## Installation

The installation has been tested on Debian and Ubuntu Linux.

1. Install dependency packages:

```
sudo apt install build-essential git zlib1g-dev
```

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

### Compile on Mac OS


1. Compile

```
cd cumulus_feature_barcoding
make all
```

## Usage
