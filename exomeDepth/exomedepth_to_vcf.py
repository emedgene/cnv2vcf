import argparse
import csv
import sys

from pysam import FastaFile

FORMAT = 'DP:BF:RR'
INFO_INDEX = 7
CNV_TYPE_TO_SHORT = {
    'deletion': 'DEL',
    'duplication': 'DUP'
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-path', dest='input_path', required=True, help='Path for input file.')
    parser.add_argument('-o', '--output-path', dest='output_path', required=True, help='Path for output file.')
    parser.add_argument('-r', '--genome-ref', dest='genome_ref', required=False,
                        default='/opt/data/ref/Human/Hg19/genome_ref/hg19.fa', help='Path to genome reference.')
    return parser.parse_args()


def get_vcf_headers(sample_name, genome_ref):
    cmdline = " ".join(sys.argv)
    headers = ['##fileformat=VCFv4.1', '##source=exomeDepthVCFConverter,version=0.1.0',
               f'##ConverterCMDLine={cmdline}']
    headers += [f"##contig=<ID={name},length={length}>"
                for name, length in zip(genome_ref.references, genome_ref.lengths)]
    headers.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}")
    return headers


def extract_sample_name(bed_file_path):
    with open(bed_file_path) as bed_file:
        temp_reader = csv.DictReader(bed_file, delimiter='\t')
        cnv_line = next(temp_reader)
        return cnv_line['sample']


def main(args):
    sample_name = extract_sample_name(args.input_path)
    with open(args.input_path) as cnv_input, FastaFile(args.genome_ref) as genome_ref,\
            open(args.output_path, 'w') as vcf_output:
        is_full_chrom_name = genome_ref.references[0].startswith('chr')
        cnv_reader = csv.DictReader(cnv_input, delimiter='\t')
        vcf_output.write('\n'.join(get_vcf_headers(sample_name, genome_ref)) + '\n')
        for cnv_line in cnv_reader:
            vcf_line = get_vcf_line(cnv_line, genome_ref, is_full_chrom_name)
            vcf_output.write(vcf_line + '\n')


def info_dict_to_string(info_dict):
    return ';'.join(sorted(("%s=%s" % info_pair for info_pair in info_dict.items()))).replace(" ", "_")


def get_sv_type(cnv_line):
    cnv_type = cnv_line['type']
    short_type = CNV_TYPE_TO_SHORT[cnv_type]
    return short_type


def get_cnv_info(cnv_line):
    num_calls = cnv_line["num.calls"]
    reads_expected = cnv_line["reads.expected"]
    start_p = cnv_line["start.p"]
    end_p = cnv_line["end.p"]
    num_exons = cnv_line["nexons"]
    sv_type = get_sv_type(cnv_line)
    start_position = cnv_line["start"]
    end_position = cnv_line["end"]
    sv_len = int(end_position) - int(start_position)
    cnv_info = {
        'SVLEN': sv_len,
        'END': end_position,
        'SVTYPE': sv_type,
        'num_calls': num_calls,
        'reads_expected': reads_expected,
        'start_p': start_p,
        'end_p': end_p,
        'num_exons': num_exons
    }
    return cnv_info


def get_sample_data(cnv_line):
    depth = cnv_line['reads.observed']
    bf = cnv_line['BF']
    read_ratio = cnv_line['reads.ratio']
    return f'{depth}:{bf}:{read_ratio}'


def get_vcf_line(cnv_line, genome_ref, is_full_chrom_name):
    chrom = get_chrom(cnv_line)
    start_position = cnv_line["start"]
    chrom_for_ref = chrom
    if is_full_chrom_name and not chrom.startswith('chr'):
        chrom_for_ref = f'chr{chrom_for_ref}'
    if not is_full_chrom_name and chrom.startswith('chr'):
        chrom_for_ref = chrom_for_ref[3:]
    ref_allele = genome_ref.fetch(region='{chr}:{pos}:{pos}'.format(chr=chrom_for_ref, pos=start_position))
    alt = get_alt(cnv_line)
    info = get_cnv_info(cnv_line)
    sample_data = get_sample_data(cnv_line)
    vcf_fields = [chrom, start_position, '.', ref_allele.upper(), alt, '.', 'PASS', info, FORMAT, sample_data]
    vcf_fields[INFO_INDEX] = info_dict_to_string(info)
    return "\t".join(vcf_fields)


def get_chrom(cnv_line):
    return cnv_line["chromosome"]


def get_alt(cnv_line):
    sv_type = get_sv_type(cnv_line)
    return f'<{sv_type}>'


if __name__ == '__main__':
    run_args = parse_args()
    main(run_args)
