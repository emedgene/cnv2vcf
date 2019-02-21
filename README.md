# BGformatConverter

Convert BG format to VCF v4.2 format Welcome to BGformatVCF converter! This tool is designed to convert Baylor Genetics CNV format to proper VCF v4.3 format and merge it with VCF file containing SNPs Dependencies:

    BCFTools: https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 can be either installed using command: sudo apt install bcftools or download from the link provided above and follow the instruction present in the README.md file
    tabix: https://excellmedia.dl.sourceforge.net/project/samtools/tabix/tabix-0.2.6.tar.bz2 can be either installed using command: sudo install tabix or donwload from the link provided anove and follow the instruction present in the README.md file
    argparse: this is python3 dependency and it comes automatically with python3 installation, hence not required to be installed seperately

Command required to run the tool:

usage: BGformat2VCF_v2.py [-h] -icnv INPUTCNV -o OUTPUT -s SAMPLEID -isnv INPUTSNV

optional arguments:

  -h, --help            show this help message and exit

Required arguments:

-icnv Input BG CNV file

-o output file containing both CNVs and SNVs in compressed VCF format

-s Sample ID to be searched in VCF file containing SNVs

-isnv Input VCF file containg SNVs

Input Requirement:

    input file containing CNVs: a) file should contain 1 header line containing these headings

    b) CNV file should atleast have 4 required columns in any order seperated by tab or space:

     "chr", "start", "end" and "logRatio"

    c) All the other information present in the file will be added to the info col of the VCF file

    input file containing SNV: a) file should be in VCF format

    b) it can be either compressed or uncompressed

    c) sample ID used as input(-s) must be contained in the VCF file

    output VCF file:

    a) is compressed sorted VCF file in BgZip(.gz) format
