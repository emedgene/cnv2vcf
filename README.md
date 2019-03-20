# cnv2vcf

This tool is designed to merge tsv containing CNVs with a vcf of SNVs

### Requirements
before using this tool please install the packages in `requirements.txt`

### Usage
```
usage: bg_cnv2vcf.py [-h] -icnv ICNV -o OUTPUT_PATH -p PROBAND_ID -f
                     FASTA_PATH -isnv ISNV

optional arguments:
  -h, --help      show this help message and exit

Required arguments::
  -icnv ICNV      Input CNV file
  -o OUTPUT_PATH  output file containing both CNVs and SNVs in VCF format
  -p PROBAND_ID   Sample ID of the proband
  -f FASTA_PATH   Path to fasta file
  -isnv ISNV      Input VCF file containing SNVs
```
