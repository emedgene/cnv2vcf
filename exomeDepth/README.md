# cnv2vcf - ExomeDepth

This tool converts tsv output of ExomeDepth algorithm into a vcf.

### Requirements
before using this tool please install the packages in `requirements.txt`

### Usage
```
usage: exomedepth_to_vcf.py [-h] -i INPUT_PATH -o OUTPUT_PATH [-r GENOME_REF]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_PATH, --input-path INPUT_PATH
                        Path for input file.
  -o OUTPUT_PATH, --output-path OUTPUT_PATH
                        Path for output file.
  -r GENOME_REF, --genome-ref GENOME_REF
                        Path to genome reference.
```
