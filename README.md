# VCFs to SQLite DB
## _VCFs to Metadata VCF to SQLite DB_

## Script details

- Tested on VCF version 4.2 files 
- Python based

## Checklist / Info

- Confirm all samples or VCF files are generated from same genome version (hg18/hg19) and same BED file interms of WES/CES.
- Confirm all the packages are installed in python env already.
- Based on file size and file count, the output and resource used will be heavy.
- Stat : For a set of 29 .vcf.gz files of average 3 MB size, the final .db file is 190 MB.

## Process / Steps involved

- STEP 1 : The current python script provided here will take all the .vcf files in the current folder to prepare a metadata VCF file
- STEP 2 : The metadata VCF file will be processed and converted into a .csv file as a supportive file
- STEP 3 : The modified .csv file will be converted into a SQLite .db file with primary keys as CHROM, POS and file_id.

## Libraries used

> os , subprocess, glob, pandas, re, sys, csv, sqlite3, datetime, time

## Input/Output files

- Input - .vcf.gz files
- Output - merged.vcf.gz , merged.vcf_table_stage1.1.csv, variants.db

Below are the information about output files:

> variants.db is the SQLite file as the required output
> merged.vcf.gz is the intermediate combined .vcf file, used for any analysis
> merged.vcf_table_stage1.1.csv is the input file to prepare .db output file

## Usage 

Copy the script to the location of .VCF files and run directly.
```sh
> python3 vcfs2SQLite_v1.py
```

## Model output for example :

