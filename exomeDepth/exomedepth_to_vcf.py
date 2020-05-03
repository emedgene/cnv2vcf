import argparse
import csv

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


def get_vcf_headers():
    headers = ['##fileformat=VCFv4.1', 'source=exomeDepthVCFConverter', "##contig=<ID=1,length=249250621>",
               "##contig=<ID=2,length=243199373>", "##contig=<ID=3,length=198022430>",
               "##contig=<ID=4,length=191154276>", "##contig=<ID=5,length=180915260>",
               "##contig=<ID=6,length=171115067>", "##contig=<ID=7,length=159138663>",
               "##contig=<ID=8,length=146364022>", "##contig=<ID=9,length=141213431>",
               "##contig=<ID=10,length=135534747>", "##contig=<ID=11,length=135006516>",
               "##contig=<ID=12,length=133851895>", "##contig=<ID=13,length=115169878>",
               "##contig=<ID=14,length=107349540>", "##contig=<ID=15,length=102531392>",
               "##contig=<ID=16,length=90354753>", "##contig=<ID=17,length=81195210>",
               "##contig=<ID=18,length=78077248>", "##contig=<ID=19,length=59128983>",
               "##contig=<ID=20,length=63025520>", "##contig=<ID=21,length=48129895>",
               "##contig=<ID=22,length=51304566>", "##contig=<ID=X,length=155270560>",
               "##contig=<ID=Y,length=59373566>", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tproband"]
    return headers


def main(args):
    with open(args.input_path) as cnv_input, FastaFile(args.genome_ref) as genome_ref,\
            open(args.output_path, 'w') as vcf_output:
        vcf_output.write('\n'.join(get_vcf_headers()) + '\n')
        cnv_reader = csv.DictReader(cnv_input, delimiter='\t')
        for cnv_line in cnv_reader:
            vcf_line = get_vcf_line(cnv_line, genome_ref)
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


def get_vcf_line(cnv_line, genome_ref):
    chrom = get_chrom(cnv_line)
    start_position = cnv_line["start"]
    ref_allele = genome_ref.fetch(region='{chr}:{pos}:{pos}'.format(chr=chrom, pos=start_position))
    alt = get_alt(cnv_line)
    info = get_cnv_info(cnv_line)
    sample_data = get_sample_data(cnv_line)
    vcf_fields = [chrom, start_position, '.', ref_allele.upper(), alt, '.', 'PASS', info, FORMAT, sample_data]
    vcf_fields[INFO_INDEX] = info_dict_to_string(info)
    return "\t".join(vcf_fields)


def get_chrom(cnv_line):
    chrom = cnv_line["chromosome"]
    if not chrom.startswith('chr'):
        chrom = 'chr' + chrom
    return chrom


def get_alt(cnv_line):
    sv_type = get_sv_type(cnv_line)
    return f'<{sv_type}>'


if __name__ == '__main__':
    run_args = parse_args()
    main(run_args)
