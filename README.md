# Cumulus Tool on Feature Barcoding

[![](https://img.shields.io/github/v/release/lilab-bcb/cumulus_feature_barcoding.svg)](https://github.com/lilab-bcb/cumulus_feature_barcoding/releases)

A fast C++ tool to extract feature-count matrix from sequence reads in FASTQ files. It is used by Cumulus for feature-count matrix generation of cell hashing, nucleus hashing, CITE-Seq and Perturb-seq protocols, using either 10x Genomics V2 or V3 chemistry.

## Installation

The installation has been tested on Linux and MacOS.

1. Install dependency packages:

On MacOS:

```
brew install cmake git
```

On Debian or Ubuntu Linux:

```
sudo apt install build-essential git
```

2. Install dependency packages:

```
git clone https://github.com/lilab-bcb/cumulus_feature_barcoding.git
```

3. Enter the directory:

```
cd cumulus_feature_barcoding
```

4. Download FQFeeder:

```
git clone https://github.com/rob-p/FQFeeder.git
```

5. Compile:

```
mkdir build
cmake ..
```

6. Now you'll have an executable named ``generate_count_matrix_ADTs`` inside your build folder. Type

```
./build/generate_count_matrix_ADTs
```

to see its usage.


## Usage
