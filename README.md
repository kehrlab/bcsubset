
bcsubset - BAM file subsetting by barcode 
=========================================

This tool selects a subset of records from a BAM file based on a barcode whitelist. 

<b>Input</b>: A BAM file and a barcode whitelist.

<b>Output</b>: A BAM file containing only  the records with the barcodes in the whitelist.

## Quickstart

### Installation
```
git clone https://github.com/kehrlab/bcsubset.git
cd bcsubset/
make
```

###  Execution
Basic execution using default values for trimming (Default: 0) and barcode tag (Default: CB):
``` 
bcsubset -w myWhitelist.txt -o outBamName.bam myBam.bam
```
Setting custom values for trimming and barcode tag:
``` 
bcsubset -w myWhitelist.txt -o outBamName.bam -t myNum -b myTag myBam.bam
```
Note: The whitelist file must contain each barcode in a new line.

## Dependencies for Installation via Make

bcsubset has the following dependencies:

-   GCC 4.9 or higher
-   ZLIB

ZLIB can simply be installed using apt:
```
sudo apt install zlib1g-dev
```

## Version and License
```
    Last update: 2021-06-24
    ctProcess version: 0.0.1
    SeqAn version: 2.4.0
    Author: Ana Belen Solis Pinson (ana.pinson[at]bihealth.de)
