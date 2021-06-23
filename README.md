bcsubset - BAM file subsetting by barcode 
=========================================

This tool selects a subset of records from a BAM file based on a barcode whitelist

<b>Input</b>: A BAM file and a barcode whitelist

<b>Output</b>: A BAM file containing only  the records with the barcodes in the whitelist.

## Quickstart

### Installation
```
git clone https://github.com/kehrlab/bcsubset.git
cd bcsubset/
make
```

###  Execution
``` 
bcsubset myBam.bam -f myWhitelist.txt -o outBamName.bam
```
Note: The whitelist text file must contain each barcode in a new line.

## Version and License
```
    Last update: 2021-06-23
    ctProcess version: 0.0.1
    SeqAn version: 2.4.0
    Author: Ana Belen Solis Pinson (ana.pinson[at]bihealth.de)
```
