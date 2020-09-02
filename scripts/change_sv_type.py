#!/usr/bin/env python

"""
This script changes the ID of a VCF file
Author: Alba Sanchis Juan (as2635@cam.ac.uk)
"""

import argparse
import vcf
import sys
from collections import Counter
from itertools import takewhile

def get_items_upto_count(dct, n):
    data = dct.most_common()
    val = data[n-1][1] #get the value of n-1th item
    return(list(takewhile(lambda x: x[1] >= val, data)))


def main():
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(description='Assign SVtype to intra-chromosomal BND')
    parser.add_argument('--vcf', dest='vcf', help='VCF input file')
    parser.add_argument('--out', dest='output', help='Output file')
    args = parser.parse_args()
    vcfin = args.vcf
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

    # samp = str(vcf_reader.samples[0])

    for record in vcf_reader:

        vtype = str(record.INFO['SVTYPE'])
        valt = str(record.ALT)
        vchr1 = str(record.CHROM)
        vchr2 = str(record.INFO["CHR2"])

        # If intrachromosomal BND - change
        if not valt.startswith('<') and vchr1 == vchr2:
            # new_alt = ''.join(['<', vtype, '>'])
            new_alt = ['<' + vtype + '>']
            record.ALT = new_alt

        if vtype not in ['BND', 'TRA']:

            ori_ids = record.samples[0]['ID']
            sur_types = record.samples[0]['TY']

            svcallers = [e.split('_')[5] for e in ori_ids.split(';')]
            svtypes = [e.split('_')[4] for e in ori_ids.split(';')]

            sur_svtypes = sur_types.split(',')

            dtypes = dict(zip(svcallers, svtypes))

            #keep sur type for survivor cx types
            if 'sniffles' in dtypes:
                if dtypes.get('sniffles') not in ['DEL', 'DUP', 'INS', 'INV']:
                    snif_pos = list(dtypes).index('sniffles')
                    dtypes['sniffles'] = sur_svtypes[snif_pos]

            dtypes_nobnd = {key: val for key, val in dtypes.items() if val != "BND"}

            #determine highest frequency type if not only BND called
            #else leave the INV type
            if len(dtypes_nobnd) > 0:
                counts = Counter(dtypes_nobnd.values())
                counts_mult = get_items_upto_count(counts, 1)

                #if more than one caller have the max type count
                #keep by order: nanosv>manta>sniffles>svim>canvas
                if len(counts_mult) == 1:
                    record.INFO['SVTYPE'] = str(counts_mult[0][0])
                else:
                    types_max = [item[0] for item in counts_mult]
                    dtypes_types_max = {key: val for key, val in dtypes_nobnd.items() if val in types_max}

                    prior_caller = ('nanosv', 'Manta', 'sniffles', 'svim', 'Canvas')

                    ordered_dtypes_max = [(key, dtypes_types_max[key]) for key in prior_caller if key in dtypes_types_max]

                    new_type = ordered_dtypes_max[0][1]

                    record.INFO['SVTYPE'] = str(new_type)

        #Change ID for all merged IDs
        record.ID = record.samples[0]['ID']

        vcf_writer.write_record(record)


if __name__ == '__main__':
    main()
