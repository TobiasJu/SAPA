#!/usr/bin/env python
# -*- coding: utf-8 -*-

# this script takes a ClinVar tabular text file with SNPs and converts it to the needed SAPA format.

import csv
import sys
import argparse

# argparse for information
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", help="clinvar file")
parser.add_argument("-c", "--counter", type=int, help="select number of SNPs to convert (default all)")
args = parser.parse_args()

if not args.input_file:
    print "ERROR, please enter an input file"
    parser.print_help()
    sys.exit(0)

print "converting file: " + args.input_file + "  to SAPA format"

counter = 0
clinvar = open(args.input_file + '.csv', 'w')
clinvar.write('"Chr","Pos","Ref","Alt","Clinical significance","Gene"\n')
with open(args.input_file) as csvfile:
    variant_lines = csv.reader(csvfile, delimiter='\t', quotechar='"')
    # skip header if there
    has_header = csv.Sniffer().has_header(csvfile.read(100))
    csvfile.seek(0)  # rewind
    incsv = csv.reader(csvfile)
    if has_header:
        header = next(variant_lines)
        print header
    for variant_line in variant_lines:
        export_string = ""
        ref_alt = variant_line[0].split(">")
        if len(ref_alt) == 2:
            ref = ref_alt[0][-1]
            alt = ref_alt[1][0]
            chr = variant_line[6]
            pos = variant_line[7]
            clin_sig = variant_line[4].split("(")
            buildver = variant_line[8]
            gene = variant_line[1]

            if buildver == "GRCh38" and chr and pos:  # and clin_sig == "Pathogenic"
                counter += 1
                print chr, pos, ref, alt, clin_sig[0], gene
                print counter
                # print buildver
                export_string += '"'+chr+'"'+","+'"'+pos+'"'+","+'"'+ref+'"'+","+'"'+alt+'"'+","+'"'+clin_sig[0]+\
                                 '"'+","+'"'+gene+'"'+"\n"
                clinvar.write(export_string)

        if args.counter:
            if counter == args.counter:
                break
clinvar.close()
print "END"