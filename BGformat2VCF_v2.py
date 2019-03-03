import sys
import os
import argparse


def write_headers(args, info_field, output_file):
    # TO ADD ALL INFORMATION FROM THE INPUT BG FILE:
    output_file.write("##fileformat=VCFv4.2\n")
    description = "Structural variant type for for this record. DUP for duplication nad DEL for deletion"
    output_file.write(f'##INFO=<ID=SVTYPE,Number=.,Type=String,Description="{description}">\n')
    output_file.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of CNV">\n')
    output_file.write(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
    output_file.write('##INFO=<ID=LogRatio,Number=1,Type=Float,Description="Log Ration for the CNV">\n')
    output_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    for info_field1 in info_field:
        line_to_be_written = '##INFO=<ID=' + info_field1 + ',Number=.,Type=String,Description="' + info_field1 + '">' + "\n"
        output_file.write(line_to_be_written)
    output_file.write('##ALT=<ID=NOCNV,Description="NO CNV found as per log-ratio(i.e 0)">\n')
    output_file.write('##ALT=<ID=DEL,Description="Deletion">\n')
    output_file.write('##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">\n')
    output_file.write('##ALT=<ID=DEL:ME:L1,Description="Deletion of L1 element">\n')
    output_file.write('##ALT=<ID=DUP,Description="Duplication">\n')
    output_file.write('##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">\n')
    output_file.write('##ALT=<ID=INS,Description="Insertion of novel sequence">\n')
    output_file.write('##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">\n')
    output_file.write('##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">\n')
    output_file.write('##ALT=<ID=INV,Description="Inversion">\n')
    output_file.write('##contig=<ID=1,length=249250621>\n')
    output_file.write('##contig=<ID=2,length=243199373>\n')
    output_file.write('##contig=<ID=3,length=198022430>\n')
    output_file.write('##contig=<ID=4,length=191154276>\n')
    output_file.write('##contig=<ID=5,length=180915260>\n')
    output_file.write('##contig=<ID=6,length=171115067>\n')
    output_file.write('##contig=<ID=7,length=159138663>\n')
    output_file.write('##contig=<ID=8,length=146364022>\n')
    output_file.write('##contig=<ID=9,length=141213431>\n')
    output_file.write('##contig=<ID=10,length=135534747>\n')
    output_file.write('##contig=<ID=11,length=135006516>\n')
    output_file.write('##contig=<ID=12,length=133851895>\n')
    output_file.write('##contig=<ID=13,length=115169878>\n')
    output_file.write('##contig=<ID=14,length=107349540>\n')
    output_file.write('##contig=<ID=15,length=102531392>\n')
    output_file.write('##contig=<ID=16,length=90354753>\n')
    output_file.write('##contig=<ID=17,length=81195210>\n')
    output_file.write('##contig=<ID=18,length=78077248>\n')
    output_file.write('##contig=<ID=19,length=59128983>\n')
    output_file.write('##contig=<ID=20,length=63025520>\n')
    output_file.write('##contig=<ID=21,length=48129895>\n')
    output_file.write('##contig=<ID=22,length=51304566>\n')
    output_file.write('##contig=<ID=X,length=155270560>\n')
    output_file.write('##contig=<ID=Y,length=59373566>\n')
    output_file.write('##contig=<ID=MT,length=16569>\n')
    line_to_be_written = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + args.s + "\n"
    output_file.write(line_to_be_written)


def parse_arguments():
    argument_parser = argparse.ArgumentParser()
    required_named = argument_parser.add_argument_group('Required arguments:')
    required_named.add_argument("-icnv", help="Input BG CNV file", required=True)
    required_named.add_argument("-o", help="output file containing both CNVs and SNVs in compressed VCF format",
                                required=True)
    required_named.add_argument("-s", help="Sample ID to be searched in VCF file containing SNVs", required=True)
    required_named.add_argument("-isnv", help="Input VCF file containg SNVs", required=True)
    args = argument_parser.parse_args()
    return args


def main():
    args = parse_arguments()

    input_file = open(str(args.icnv), 'r')
    output_file = open(".temp.vcf", 'w')

    # TO REMOVE '#' FROM HEADER LINE(IF PRESENT)
    header_line = input_file.readline()
    if header_line.startswith('#'):
        header_line = header_line.lstrip('#')
    headers = header_line.lower().split()

    info_field = list(set(headers) - {'chr', 'seqnames', 'start', 'end', 'logratio'})
    info_field_index = {header: index for index, header in enumerate(headers)}

    info_field_index = []
    # TO EXTRACT ALL HEADER INFORMATION FROM THE INPUT BG FILE
    for index, header_variable in enumerate(headers):
        if "chr" == header_variable:
            chr_index = index
        elif "seqnames" == header_variable:
            chr_index = index
        elif "start" == header_variable:
            start_index = index
        elif "end" == header_variable:
            end_index = index
        elif "logratio" == header_variable:
            logratio_index = index
        else:
            info_field_index.append(index)

    write_headers(args, info_field, output_file)

    for line in input_file.readlines():
        line_elements = line.split()
        # TO CHECK TYPE OF CNV:
        if float(line_elements[logratio_index]) == 0.0:
            cnv_type = '<NOCNV>'
            genotype == '0/0'
        elif float(line_elements[logratio_index]) > 0.0:
            cnv_type = '<INS>'
            if 0 < float(line_elements[logratio_index]) < 0.5:
                genotype = '0/1'
            else:
                genotype = '1/1'
        else:
            cnv_type = '<DEL>'
            if 0 > float(line_elements[logratio_index]) > -0.5:
                genotype = '0/1'
            elif -0.51 > float(line_elements[logratio_index]):
                genotype = '1/1'

        cnv_length = int(line_elements[end_index]) - int(line_elements[start_index]) + 1

        # TO CHECK IF END COORDIANTE IS LESS THAN STAR COORDINATE:
        if cnv_length < 0:
            print("\nERROR: End coordiante is less than Start coordiante for chr", line_elements[chr_index], ":",
                  line_elements[start_index], "-", line_elements[end_index], sep="")
            sys.exit(1)

        # ERROR HANDLING IN CASE EITHER OF CHR or START or END or LOGRATIO NOT PRESENT IN THE BG FILE
        try:
            logratio_index
            chr_index
            start_index
            end_index
        except NameError:
            logration_index = None
            chr_index = None
            start_index = None
            end_index = None
        if logratio_index is None:
            print("logRatio not defined in the input file header")
            sys.exit(0)
        elif chr_index is None:
            print("chr not defined in the input file header")
            sys.exit(0)
        elif start_index is None:
            print("start not defined in the input file header")
            sys.exit(0)
        elif end_index is None:
            print("end not defined in the input file header")
            sys.exit(0)
        else:
            #		if not line_elements[chr_index].startswith("chr"):  #ONLY REQUIRED IF USER WANT TO ADD CHR AT THE START OF EACH CHROMOSOME
            #			output.write("chr",end="")
            line_to_be_written = line_elements[chr_index] + "\t" + str(line_elements[start_index]) + "\t.\t.\t" + str(
                cnv_type) + "\t.\t.\tSVTYPE=" + str(cnv_type) + ";SVLEN=" + str(cnv_length) + ";END=" + str(
                line_elements[end_index]) + ";LogRatio=" + str(line_elements[logratio_index])
            output_file.write(line_to_be_written)

        # TO ADD ALL EXTRA INFORMATION TO INFO FIELD OF VCF
        i = 0
        for info_field_var in info_field:
            i1 = info_field_index[i]
            line_to_be_written = ";" + info_field_var + "=" + line_elements[i1]
            output_file.write(line_to_be_written)
            i += 1
        line_to_be_written = "\tGT\t" + genotype + "\n"
        output_file.write(line_to_be_written)

    input_file.close()
    output_file.close()

    # TO SORT CNV OUTPUT IN VCF FORMAT:
    line_to_be_executed = "bcftools sort -O v -o .temp1.vcf .temp.vcf"
    os.system(line_to_be_executed)

    # TO COMPRESS CNV VCF FILE IN BGZIP FORMAT:
    line_to_be_executed = "bgzip -c .temp1.vcf > .temp1.vcf.gz"
    os.system(line_to_be_executed)

    # TO TABIX INDEX CNV VCF FILE:
    line_to_be_executed = "tabix -p vcf .temp1.vcf.gz"
    os.system(line_to_be_executed)

    # TO EXTRACT SAMPLE FROM SNV VCF FILE:
    line_to_be_executed = "bcftools view -s " + args.s + " -o .snv_temp.vcf.gz -O z " + args.isnv
    os.system(line_to_be_executed)

    # TO TABIX INDEX SNV VCF FILE:
    line_to_be_executed = "tabix -p vcf .snv_temp.vcf.gz"
    os.system(line_to_be_executed)
    # TO CONCATNATE CNV+SNV VCF FILES TOGETHER:
    line_to_be_executed = "bcftools concat -a -o " + args.o + " -O z .snv_temp.vcf.gz .temp1.vcf.gz"
    os.system(line_to_be_executed)

    print("\n\tFINAL OUTPUT compressed (SNP+CNV) VCF FILE: ", args.o, "\n\n")
    # os.system("rm .temp.vcf .temp1.vcf .temp1.vcf.gz .snv_temp.vcf.gz")


if __name__ == '__main__':
    main()
