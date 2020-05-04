#!/usr/bin/env python

"""
This script changes the ID of a VCF file
Author: Alba Sanchis Juan (as2635@cam.ac.uk)
"""

import argparse
import re
import vcf
import sys


def main():
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(description='Assign SVtype to intra-chromosomal BND')
    parser.add_argument('--vcf', dest='vcf', help='VCF input file')
    parser.add_argument('--caller', dest='caller', help='Variant caller algorithm')
    parser.add_argument('--out', dest='output', help='Output file')
    args = parser.parse_args()
    vcfin = args.vcf
    caller = args.caller
    output = args.output

    """
    Read files
    """

    if vcfin is '-':
        vcf_reader = vcf.Reader(sys.stdin)
    else:
        vcf_reader = vcf.Reader(open(vcfin, 'r'))

    if output is None:
        vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
    else:
        vcf_writer = vcf.Writer(open(output, 'w'), vcf_reader)

    """
    Edit ID field
    """
    for record in vcf_reader:

        vtype = str(record.INFO["SVTYPE"])
        vchr1 = str(record.CHROM)
        vstart = str(record.POS)
        valt = str(record.ALT)

        if vtype == "DUP:TANDEM":
            record.INFO["SVTYPE"] = "DUP"
            vtype = "DUP"

        if vtype == "BND":
            coords = re.split('\]|\[', valt)[2]
            vchr2 = re.split(':', coords)[0]
            vend = re.split(':', coords)[1]
        else:
            vchr2 = str(record.CHROM)
            vend = str(record.INFO["END"])

        new_id = '_'.join([vchr1, vstart, vchr2, vend, vtype, caller])

        record.ID = new_id

        if vstart != "0":
            if vchr1 == vchr2:
                vlen = (int(vend) - int(vstart)) + 1
                if vlen >= 50:
                    vcf_writer.write_record(record)
            else:
                vcf_writer.write_record(record)


if __name__ == '__main__':
    main()
