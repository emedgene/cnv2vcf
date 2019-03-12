import argparse
import csv
import gzip
from contextlib import contextmanager

from Bio import bgzf
import sys
from pysam import FastaFile


SAMPLE_DATA_FIELDS = ['GT', 'ZS', 'LR', 'CD', 'RPKM']


@contextmanager
def open_file(filename, mode='r'):
    if 'w' in mode and filename.endswith('.gz'):
        with bgzf.open(filename, mode) as f:
            yield f
    else:
        if 't' not in mode:
            mode += 't'
        _open = gzip.open if filename.endswith('.gz') else open
        with _open(filename, mode, encoding='utf-8') as f:
            yield f


def parse_arguments():
    argument_parser = argparse.ArgumentParser()
    parser = argument_parser.add_argument_group('Required arguments:')
    parser.add_argument("-icnv", help="Input BG CNV file", required=True)
    parser.add_argument("-o", dest='output_path', required=True,
                        help="output file containing both CNVs and SNVs in compressed VCF format")
    parser.add_argument("-p", dest='proband_id', help="Sample ID of the proband", required=True)
    parser.add_argument("-f", help="Path to fastq file", dest='fasta_path', required=True)
    parser.add_argument("-isnv", help="Input VCF file containg SNVs", required=True)
    return argument_parser.parse_args()


def get_cnv_type(log_ratio):
    return 'INS' if log_ratio > 0 else 'DEL' if log_ratio < 0 else 'NOCNV'


def main():
    args = parse_arguments()
    with open_file(args.output_path) as output_file, open_file(args.icnv) as input_cnv_file,\
            FastaFile(args.fasta_path) as genome_ref, open_file(args.isnv) as vcf_input:
        # Write all snv lines
        num_samples = 0
        proband_index = 0
        for snv_line in vcf_input:
            if snv_line.startswith('#'):
                if not snv_line.startswith('##'):
                    split_line = snv_line.split("\t")
                    num_samples = len(split_line) - 8  # chr, pos, id, ref, alt, qual, filter, info, format
                    try:
                        proband_index = split_line.index(args.proband_id)
                    except ValueError:
                        print("Sample ID of proband does not exist in VCF header")
                        sys.exit(2)
            output_file.write(snv_line)
        if not num_samples or not proband_index:
            print("Missing VCF header line")
            sys.exit(3)

        cnv_reader = csv.DictReader(input_cnv_file, delimiter='\t')
        for cnv_line in cnv_reader:
            log_ratio = float(cnv_line['LogRatio'])
            cnv_type = get_cnv_type(log_ratio)
            cnv_length = int(cnv_line['end']) - int(cnv_line['start']) + 1

            chrom = cnv_line['seqnames']
            position = cnv_line['start']
            ref_allele = genome_ref.fetch(region='{chr}:{pos}:{pos}'.format(chr=chrom, pos=position))
            info_dict = {'SVTYPE': cnv_type, 'SVLEN': cnv_length, 'END': cnv_line['end'], 'LogRatio': log_ratio}
            info_field = ';'.join([f'{key}={value}' for key, value in info_dict.items()])

            # chr, pos, id, ref, alt, qual, filter, info, format, *samples
            vcf_fields = [chrom, position, '.', ref_allele, f'<{cnv_type}>', '.', '.', info_field,
                          ':'.join(SAMPLE_DATA_FIELDS)]
            vcf_fields += ['./.' + ':'.join(['.'] * (len(SAMPLE_DATA_FIELDS) - 1))] * num_samples
            vcf_fields[proband_index] = ':'.join(['1/1', cnv_line['Zscore'], cnv_line['LogRatio'],
                                                  cnv_line['CountDensity'], cnv_line['RPKM']])
            output_file.write("\t".join(vcf_fields) + "\n")


if __name__ == '__main__':
    main()
