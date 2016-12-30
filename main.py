#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

try:
    import setuptools
except ImportError, e:
    print "PiP is not installed"
    sys.exit(1)

try:
    import dominate
except ImportError, e:
    pass
    print "Dominate is missing"
    sys.exit(1)

from dominate.tags import *
from os.path import exists
import re
import csv
import os
import collections
import subprocess
import ensembl_rest
import itertools
from snp import SNP
from allProt import AllProt
from probedMutation import ProbedMutation
from dna_to_aa_translator import Translator
from geneDNA import GeneDNA
from annovarParser import AnnovarParser
import argparse

# argparse for information
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
# group.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
group.add_argument("-q", "--quiet", action="store_true", help="prevent output in command line")
# parser.add_argument("-l", "--log", action="store_true", help="store the output in a log file")
# parser.add_argument("-n", "--number", type=int, help="just a Test number")
parser.add_argument("-m", "--manual", action="store_true", help="display the manual for this program")
parser.add_argument("-i", "--input_file", help="comma separated table with SNP's")
parser.add_argument("-d", "--detail", action="store_true", help="write detailed output file")
parser.add_argument("-s", "--separator", help='set the input file separator (default: ",")')
parser.add_argument("-t", "--text_delimiter", help='set the input text delimiter (default: ")')
parser.add_argument("-b", "---buildver", type=str, help="database buildversion (default: hg19)")
parser.add_argument("-o", "--output_file", type=str, help="output file name (default output.csv)")
parser.add_argument("-f", "--fast", action="store_true", help="run annotation just with a region based approach, "
                                                              "for faster computing and less download file demand")
parser.add_argument("-fi", "--filter", action="store_true", help="filter SNPs for nonsynonymus and clinically "
                                                                 "significant (>95%%) SNPs")
parser.add_argument("-fd", "--filter_deleterious", action="store_true", help="output only deleterious SNPs")
args = parser.parse_args()

# sanity check ###
if not len(sys.argv) > 1:
    parser.print_help()
    sys.exit(0)

if args.manual:
    print """SAPA - SNP annotation programm for AML

75% of a Acute myeloid leukemia (AML) patients SNPs are unique to them and their impact is still unresolved.

SAPA is designed to add additional information to an Illumina truseq amplicon variants csv file. Such as Scores
for nonsynonymous scoring matrices, divided into 3 different combined Scores, function prediction scores,
conservation scores and one ensemble score. So that the viewer gets more Information whether a mutation is delirious
or harmless. Due to performance and data reasons, SAPA uses annovar1 and their precomputed SNP database, to quickly
gather the deleterious prediction methods scores.

First run:

./main.py -i data/truseq_example_data.csv

Please Note that during the first run of the program, the required database will be downloaded. This may take some time,
 depending on your internet connection.

Parameters:

-h, --help : show the help message and exit

-q, --quiet : prevent output in command line, run the program and don't bother

-m, --manual : display the manual for this program

-i, --input_file : truseq amplicon table containing SNPs
    > e.g. -i data/truseq_example_data.csv

-d, --detail : write detailed output file, with all the single scores of the deleteriousness prediction methods
               for nonsynonymous SNVs and amino acid substitution (if any)

-s, --separator : set the input file separator (default: ",")
    > e.g. -s ";" this will set the column separator to a semicolon

-t, --text_delimiter : set the input text delimiter (default: ")
    > e.g. -t "'" this will set the text delimiter to an apostrophe, if a column contains multiple entries

-id, --input_directory : hg19 (GRCh37) database directory (default /hg19)
    > e.g. -id humandb/ set the input directory to be named humandb

-o, --output_file : output file name (default output.csv)
    > e.g. -o annotated_snps_detailed.txt -d running the detailed operation mode, and saving it in the file
              annotated_snps_detailed.txt

-f, --fast : run annotation just with a region based approach, for faster computing and less download file demand

-fi, --filter : filter all entries and only use nonsynonymus and clinically significant (>95%) SNPs in the output file

Several commonly used databases are integrated: ‘cytoBand’ for the chromosome coordinate of each cytogenetic band,
‘exac03’ for the variants reported in the Exome Aggregation Consortium (version 0.3)50, ‘dbnsfp30a’ for various
functional deleteriousness prediction scores from the dbNSFP database (version 2.6)51,
‘clinvar_20140929’ for the variants reported in the ClinVar database (version 20140929)52 and
‘avsnp147’ for the dbSNP database (version 138).
             """
    parser.print_help()
    sys.exit(0)

if not args.input_file:
    print "ERROR, please enter a input file parameter"
    parser.print_help()
    sys.exit(0)

if not args.quiet:
    print args

if args.buildver == "hg38" or args.buildver == "GRCh38":
    buildversion = "hg38"
    print "buildversion is GRCh38/hg38"
elif args.buildver == "dm3" or args.buildver == "mm9":
    print "Please contact me via github and tell my why, if you need different build versions."
else:
    buildversion = "hg19"
    print "buildversion is GRCh37/hg19"

# load the input file into a mutation Objects from the given patient Data ###
fail_count = 0
success_count = 0
with open(args.input_file) as csvfile:
    if args.separator and args.text_delimiter:
        variant_lines = csv.reader(csvfile, delimiter="{}".format(args.separator),
                                   quotechar="{}".format(args.text_delimiter))
    elif args.separator and not args.text_delimiter or args.text_delimiter and not args.separator:
        print "please enter delimiter AND quote chars"
        sys.exit(0)
    else:
        variant_lines = csv.reader(csvfile, delimiter=',', quotechar='"')

    # skip header if there
    has_header = csv.Sniffer().has_header(csvfile.read(25))  # 100
    csvfile.seek(0)  # rewind
    incsv = csv.reader(csvfile)
    if has_header:
        next(variant_lines)
        if not args.quiet:
            print "skipping header"
    l_count = 0
    snps = []
    for variant_line in variant_lines:
        # remove /n form end of line
        # variant_line = variant_line.strip()
        # remove " from lines
        # variant_line[0] = variant_line[0].translate(None, '"')
        # variant_line[-1] = variant_line[-1].translate(None, '"')
        if len(variant_line) >= 16:
            success_count += 1
            context = variant_line[5].split(",")
            consequences = variant_line[6].split(",")
            snps.append(SNP(l_count, variant_line[0], int(variant_line[1]), variant_line[2], variant_line[3],
                            variant_line[4], context, consequences, variant_line[7], variant_line[8], variant_line[9],
                            int(variant_line[10]), variant_line[11], variant_line[12], variant_line[13],
                            variant_line[14], variant_line[15]))
        elif 16 > len(variant_line) > 4:
            success_count += 1
            # context = variant_line[5].split(",")
            # consequences = variant_line[6].split(",")
            snps.append(SNP(l_count, variant_line[0], int(variant_line[1]), variant_line[2], variant_line[3],
                            variant_line[4], variant_line[5], ".", ".", ".", ".", ".", ".", ".", ".", ".", "."))
        else:
            print "INVALID DATA (16 >= length > 4) in Line {}".format(l_count)
            print variant_line
            fail_count += 1
        l_count += 1
if not args.quiet:
    print "created patient SNP objects with " + str(len(snps)) + " unique SNPs\n"
csvfile.close()

if fail_count > success_count:
    print "INPUT ERROR! wrong separator selected"
    sys.exit(0)

# filter all entries in the patient mutation data set ###
if args.filter:
    print "pre filter SNP count: " + str(len(snps))
    i = 0
    for snp in snps:
        # only clinically relevant quality
        if snp.get_qual() <= 95:
            del snps[i]

        # if mutation does not change the amino acid, it does not affect the cell (in most cases)
        if "synonymous_variant" in snp.get_consequences():
            del snps[i]
        i += 1
    print "past filter SNP count: " + str(len(snps))

# write tab delimited file for annovar #
tab_mutations = open('amplicon_variants_tab.csv', 'w')
for snp in snps:
    # check if type is deletion, correction of the data for annovar
    if "Deletion" in snp.get_type():
        # print snp.get_ref()
        # print mutation.get_alt()
        newEnd = snp.get_pos() + (len(snp.get_ref()) - 2)
        snp.set_alt("-")
        if len(snp.get_ref()) > 2:
            snp.set_ref(0)
        else:
            try:
                snp.set_ref(snp.get_ref()[1])
            except IndexError:
                snp.set_ref(snp.get_ref()[0])
        snp.set_new_end(newEnd)
        # print mutation.get_pos()
        # print newEnd
        # print mutation.get_ref()
        # print snp.get_alt()
    tab_mutations.write(snp.export())
tab_mutations.close()
if not args.quiet:
    print "created tab delimited file for annovar"

# WORKS JUST UNDER UBUNTU OR THE UBUNTU BASH FOR WINDOWS #
# get annovar databases if needed ###

annotate_variation = "./perl/annotate_variation.pl "
databases = ["-buildver " + buildversion + " -downdb -webfrom annovar refGene " + buildversion,
             "-buildver " + buildversion + " -downdb cytoBand " + buildversion,
             "-buildver " + buildversion + " -downdb -webfrom annovar esp6500siv2_all " + buildversion,
             "-buildver " + buildversion + " -downdb -webfrom annovar avsnp147 " + buildversion,  # avsnp147
             "-buildver " + buildversion + " -downdb -webfrom annovar dbnsfp30a " + buildversion
             ]

if not os.path.isfile(buildversion + "/" + buildversion + "_refGene.txt"):
    print "downloading dependencies for " + buildversion + "..."
    popen = subprocess.Popen([annotate_variation + databases[0]], shell=True)
    popen.communicate()

if not os.path.isfile(buildversion + "/" + buildversion + "_refGene.txt"):
    print "download failed! Please try again or download the file directly with the link above "
    sys.exit(0)

if not os.path.isfile(buildversion + "/" + buildversion + "_dbnsfp30a.txt"):
    print "downloading dependencies for " + buildversion + "..."
    popen = subprocess.Popen([annotate_variation + databases[4]], shell=True)
    popen.communicate()

if not os.path.isfile(buildversion + "/" + buildversion + "_dbnsfp30a.txt"):
    print "download failed! Please try again or download the file directly with the link above "
    sys.exit(0)

if not args.fast:
    if not os.path.isfile(buildversion + "/" + buildversion + "_cytoBand.txt"):
        print "downloading dependencies for " + buildversion + "..."
        popen = subprocess.Popen([annotate_variation + databases[1]], shell=True)
        popen.communicate()

    if not os.path.isfile(buildversion + "/" + buildversion + "_cytoBand.txt"):
        print "download failed! Please try again or download the file directly with the link above "
        sys.exit(0)

    if not os.path.isfile(buildversion + "/" + buildversion + "_esp6500siv2_all.txt"):
        print "downloading dependencies for " + buildversion + "..."
        popen = subprocess.Popen([annotate_variation + databases[2]], shell=True)
        popen.communicate()

    if not os.path.isfile(buildversion + "/" + buildversion + "_esp6500siv2_all.txt"):
        print "download failed! Please try again or download the file directly with the link above "
        sys.exit(0)

    if not os.path.isfile(buildversion + "/" + buildversion + "_avsnp147.txt"):
        print "downloading dependencies for " + buildversion + "..."
        popen = subprocess.Popen([annotate_variation + databases[3]], shell=True)
        popen.communicate()

    if not os.path.isfile(buildversion + "/" + buildversion + "_avsnp147.txt"):
        print "download failed! Please try again or download the file directly with the link above "
        sys.exit(0)

# annotate the SNPs ###
print "annotating the SNPs"
annovar_pl = "./perl/table_annovar.pl "
dir_path = os.path.dirname(os.path.realpath(__file__))
# print dir_path

if args.fast:
    params = "amplicon_variants_tab.csv " + buildversion + " -buildver " + buildversion + " -out myanno -remove -protocol " \
                                                                                          "refGene,dbnsfp30a -operation g,f -nastring ."
else:
    params = "amplicon_variants_tab.csv " + buildversion + " -buildver " + buildversion + " -out myanno -remove -protocol " \
                                                                                          "refGene,cytoBand,esp6500siv2_all,avsnp147,dbnsfp30a " \
                                                                                          "-operation g,r,f,f,f -nastring . "  # avsnp147
if args.quiet:
    FNULL = open(os.devnull, 'w')
    popen = subprocess.Popen([annovar_pl + params], shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
else:
    print annovar_pl + params
    popen = subprocess.Popen([annovar_pl + params], shell=True)
# wait until it's finished
popen.communicate()

# parse annovar file ###
annovar = []
l_count = 0
failed = 1
a_file = "myanno." + buildversion + "_multianno.txt"  # 'myanno.hg19_multianno.txt'

with open(a_file, 'r') as annovar_file:
    for row in annovar_file:
        # filter header
        if l_count != 0:
            row = row.strip()
            row = row.split("\t")
            if len(row) == 44:  # fast run
                annovar.append(AnnovarParser(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8],
                                             row[9], ".", ".", ",", row[10], row[11], row[12], row[13], row[14],
                                             row[15], row[16],
                                             row[17], row[18], row[19], row[20], row[21], row[22], row[23], row[24],
                                             row[25], row[26], row[27], row[28], row[29], row[30], row[31], row[32],
                                             row[33], row[34], row[35], row[36], row[37], row[38], row[39], row[40],
                                             row[41], row[42], row[43]))
            elif len(row) == 47:
                annovar.append(AnnovarParser(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8],
                                             row[9], row[10], row[11], row[12], row[13], row[14], row[15], row[16],
                                             row[17], row[18], row[19], row[20], row[21], row[22], row[23], row[24],
                                             row[25], row[26], row[27], row[28], row[29], row[30], row[31], row[32],
                                             row[33], row[34], row[35], row[36], row[37], row[38], row[39], row[40],
                                             row[41], row[42], row[43], row[44], row[45], row[46]))
            else:  # fatal error
                annovar.append(AnnovarParser(".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".",
                                             ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".",
                                             ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".",
                                             ".", ".", ".", ".", ".", ".", ".", "."))
                print "ERROR could not handle this row: "
                print row
                print len(row)
                failed += 1
        l_count += 1
annovar_file.close()

if failed == l_count:
    print "FATAL ERROR: could not handle annovar input, SNP annotation failed"
    sys.exit(0)

if not args.quiet:
    print "annotated SNPs count: " + str(len(annovar))

# link annotations and SNPs together, with proofing! ###
counters = 0
snps_with_annotation = {}
for snp, annotation in zip(snps, annovar):
    if annotation._AnnovarParser__Chr == snp.get_chr() and int(annotation._AnnovarParser__Start) == snp.get_pos():
        snps_with_annotation[snp] = annotation
        counters += 1

        # iterate over annovar data and get final scores ###
        # for data in annovar:
        # score = data._AnnovarParser__SIFT_score
        # if score != ".":
        # rel_score = float(data._AnnovarParser__SIFT_score) / data._AnnovarParser__SIFT_max
        # print rel_score
###### some work to do here ####
# sys.exit(0)

# create Objects containing all human proteins ###
allHumanProteins = []
allProtFile = open('data/allprots.csv', 'r')
lines = allProtFile.readlines()[1:]
line_count = 0
for snp_entry in lines:
    # remove /n form end of line
    snp_entry = snp_entry.strip()
    variant_line = snp_entry.split('\t')
    line_count += 1
    if len(variant_line) >= 20:
        i = 0
        for split in variant_line:
            variant_line[i] = variant_line[i].translate(None, "\'")
            i += 1
        prot = variant_line[0]
        geneSyn = variant_line[1].split(",")
        ensembl = variant_line[2]
        position = variant_line[5].split("-")
        start = position[0]
        end = position[1]
        geneDesc = variant_line[3].split(",")
        chromosome = "chr" + str(variant_line[4])
        # print chromosome
        allHumanProteins.append(AllProt(prot, geneSyn, ensembl, geneDesc, chromosome, int(start),
                                        int(end), variant_line[6], variant_line[7:-1]))
allProtFile.close()
if not args.quiet:
    print "created all human protein objects"

# search after SNP corresponding genes ###
coding_mutations = []
for snp in snps:
    if "Coding" in snp.get_consequences():
        if not args.quiet:
            print "{} ID: {}, Position: {}".format("unknown mutation", snp.get_id(), snp.get_pos())
        for prot in allHumanProteins:
            if prot.get_start() < snp.get_pos() < prot.get_end() and snp.get_chr() == \
                    prot.get_chromosome():
                # print "SNP Pos: {}, Ref Gene Start: {},  Ref Gene End: {}".format(mutation.get_pos(),
                #                                                                        prot.get_start(),
                #                                                                        prot.get_end())
                if not args.quiet:
                    print "found Gene: " + prot.get_gene()
                coding_mutations.append(ProbedMutation(snp.get_id(), snp.get_chr(), snp.get_pos(),
                                                       snp.get_ref(), snp.get_alt(), snp.get_type(),
                                                       snp.get_context(), snp.get_consequences(),
                                                       snp.get_dbSNP(), snp.get_cosmic(),
                                                       snp.get_clinVar(), snp.get_qual(),
                                                       snp.get_altFreq(), snp.get_totalDepth(),
                                                       snp.get_refDepth(), snp.get_altDepth(),
                                                       snp.get_strandBias(), str(prot.get_chromosome()),
                                                       prot.get_gene(), prot.get_geneSyn(), prot.get_geneDesc(),
                                                       prot.get_proteinClass(), prot.get_start(), prot.get_end()
                                                       ))
    else:
        coding_mutations.append(ProbedMutation(snp.get_id(), snp.get_chr(), snp.get_pos(),
                                               snp.get_ref(), snp.get_alt(), snp.get_type(),
                                               snp.get_context(), snp.get_consequences(),
                                               snp.get_dbSNP(), snp.get_cosmic(),
                                               snp.get_clinVar(), snp.get_qual(),
                                               snp.get_altFreq(), snp.get_totalDepth(),
                                               snp.get_refDepth(), snp.get_altDepth(),
                                               snp.get_strandBias(), ".", ".", ".", ".", ".", ".", "."
                                               ))

# find DNA sequence for gene in each region and translate it ###
# mutation_with_sequence = {}
# for c_muta in coding_mutations:
#     print "Expected gene length: " + str(c_muta.get_geneEnd() - c_muta.get_geneStart())
#     openString = args.input_directory + "/chromFa/" + c_muta.get_geneChromosome() + ".fa"
#     hg19_chromosome = open(openString, "r")
#     with open(openString) as gf:
#         chromosome = gf.read()
#         chromosome = chromosome.replace(">" + c_muta.get_geneChromosome(), '')
#         chromosome = chromosome.replace("\n", '').replace("\r", '').replace("\n\r", '')
#
#     print len(chromosome)
#     gene = chromosome[c_muta.get_geneStart():c_muta.get_geneEnd()]
#     print len(gene)
#
#     # Translate DNA to AA #
#     amino_seq = Translator().translate_dna_sequence(gene)
#     # print amino_seq
#     g_dna = GeneDNA(c_muta.get_gene(), c_muta.get_geneChromosome(), c_muta.get_geneStart(), c_muta.get_geneEnd(),
#                     gene, amino_seq)
#
#     # print "Element name: " + element.get_name()
#
#     mutation_with_sequence[c_muta] = g_dna
#
# for c_muta, g_dna in mutation_with_sequence.iteritems():
#     print c_muta.get_gene()
#     print c_muta.get_geneChromosome()
#     print c_muta.get_pos()
#     print g_dna.get_aa_sequence()


# ensemble API for GRCh37 hg19 ###
# ensembl_rest.run(species="human", symbol="DOPEY2")


# write in export table ###
if not args.output_file:
    target = open("output.csv", 'w')
    print "writing export in: output.csv"
    target_html = open("output.html", 'w')
    print "writing html table in: output.html"
else:
    target = open(args.output_file, 'w')
    print "writing export in: " + args.output_file
    target_html = open(args.output_file + ".html", 'w')
    print "writing html table in: " + args.output_file + ".html"
export_cnt = 0
ordered_snps_with_annotation = collections.OrderedDict(sorted(snps_with_annotation.items()))

for snp, annotation in ordered_snps_with_annotation.iteritems():
    # write header first:
    if export_cnt == 0:
        if args.detail:
            header = str(snp.print_header()) + str(annotation.print_header() + "final prediction\t")
        else:
            header = str(snp.print_header()) + "Gene\tfunction prediction scores [0-1]\tconservation scores[-12.3-" \
                                               "6.17]\t" \
                                               "ensemble scores[0-60]\tfinal prediction\t"
        target.write(header)
        target.write("\n")
    # write rows in table
    try:
        altFreq = snp.get_altFreq() * 100
    except ValueError:
        altFreq = snp.get_altFreq()
    if args.detail:
        export_string = str("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
                            "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t"
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t"
                            "".format(snp.get_id(), snp.get_chr(), snp.get_pos(), snp.get_ref(), snp.get_alt(),
                                      snp.get_type(), ','.join(snp.get_context()), ','.join(snp.get_consequences()),
                                      snp.get_dbSNP(), snp.get_cosmic(), snp.get_clinVar(), snp.get_qual(),
                                      altFreq, snp.get_totalDepth(), snp.get_refDepth(),
                                      snp.get_altDepth(), snp.get_strandBias(),
                                      annotation._AnnovarParser__Chr,
                                      annotation._AnnovarParser__Start,
                                      annotation._AnnovarParser__End,
                                      annotation._AnnovarParser__Ref,
                                      annotation._AnnovarParser__Alt,
                                      annotation._AnnovarParser__Func_refGene,
                                      annotation._AnnovarParser__Gene_refGene,
                                      annotation._AnnovarParser__GeneDetail_refGene,
                                      annotation._AnnovarParser__ExonicFunc_refGene,
                                      annotation._AnnovarParser__AAChange_refGene,
                                      annotation._AnnovarParser__cytoBand,
                                      annotation._AnnovarParser__esp6500siv2_all,
                                      annotation._AnnovarParser__avsnp147,
                                      annotation._AnnovarParser__SIFT_score[0],
                                      annotation._AnnovarParser__SIFT_pred,
                                      annotation._AnnovarParser__Polyphen2_HDIV_score[0],
                                      annotation._AnnovarParser__Polyphen2_HDIV_pred,
                                      annotation._AnnovarParser__Polyphen2_HVAR_score[0],
                                      annotation._AnnovarParser__Polyphen2_HVAR_pred,
                                      annotation._AnnovarParser__LRT_score[0],
                                      annotation._AnnovarParser__LRT_pred,
                                      annotation._AnnovarParser__MutationTaster_score[0],
                                      annotation._AnnovarParser__MutationTaster_pred,
                                      annotation._AnnovarParser__MutationAssessor_score[0],
                                      annotation._AnnovarParser__MutationAssessor_pred,
                                      annotation._AnnovarParser__FATHMM_score[0],
                                      annotation._AnnovarParser__FATHMM_pred,
                                      annotation._AnnovarParser__PROVEAN_score[0],
                                      annotation._AnnovarParser__PROVEAN_pred,
                                      annotation._AnnovarParser__VEST3_score[0],
                                      annotation._AnnovarParser__CADD_raw[0],
                                      annotation._AnnovarParser__CADD_phred[0],
                                      annotation._AnnovarParser__DANN_score[0],
                                      annotation._AnnovarParser__fathmm_MKL_coding_score[0],
                                      annotation._AnnovarParser__fathmm_MKL_coding_pred,
                                      annotation._AnnovarParser__MetaSVM_score[0],
                                      annotation._AnnovarParser__MetaSVM_pred,
                                      annotation._AnnovarParser__MetaLR_score[0],
                                      annotation._AnnovarParser__MetaLR_pred,
                                      annotation._AnnovarParser__integrated_fitCons_score[0],
                                      annotation._AnnovarParser__integrated_confidence_value,
                                      annotation._AnnovarParser__GERP_RS[0],
                                      annotation._AnnovarParser__phyloP7way_vertebrate[0],
                                      annotation._AnnovarParser__phyloP20way_mammalian[0],
                                      annotation._AnnovarParser__phastCons7way_vertebrate[0],
                                      annotation._AnnovarParser__phastCons20way_mammalian[0],
                                      annotation._AnnovarParser__SiPhy_29way_logOdds[0]
                                      ))
        # snp.get_geneChromosome(),
        # snp.get_gene(), snp.get_geneSyn(), snp.get_geneDesc(),
        # snp.get_proteinClass(), snp.get_geneStart(), snp.get_geneEnd()
        # ))
    else:
        export_string = str("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t"
                            "".format(snp.get_id(), snp.get_chr(), snp.get_pos(),
                                      snp.get_ref(), snp.get_alt(), snp.get_type(),
                                      ','.join(snp.get_context()), ','.join(snp.get_consequences()),
                                      snp.get_dbSNP(), snp.get_cosmic(), snp.get_clinVar(),
                                      snp.get_qual(), snp.get_altFreq(),
                                      snp.get_totalDepth(), snp.get_refDepth(),
                                      snp.get_altDepth(), snp.get_strandBias(),
                                      annotation._AnnovarParser__Gene_refGene,
                                      annotation._AnnovarParser__MetaLR_score[0],
                                      annotation._AnnovarParser__GERP_RS[0],
                                      annotation._AnnovarParser__CADD_phred[0]
                                      ))
    # check for final prediction ###
    if annotation._AnnovarParser__MetaLR_pred == "T":
        export_string += "Tolerated"
    elif annotation._AnnovarParser__MetaLR_pred == "D":
        export_string += "Deleterious"
    elif annotation._AnnovarParser__MetaLR_pred == "P":
        export_string += "possibly damaging"
    elif annotation._AnnovarParser__MetaLR_pred == ".":
        # no Meta Socore available
        if annotation._AnnovarParser__fathmm_MKL_coding_pred == "D":
            export_string += "Deleterious (FATHMM)"
        elif not annotation._AnnovarParser__DANN_score[0] == ".":
            if float(annotation._AnnovarParser__DANN_score[0]) > 0.96:
                export_string += "Deleterious (DANN)"
        else:
            if annotation._AnnovarParser__ExonicFunc_refGene == "synonymous SNV":
                export_string += "Tolerated (synonymous)"
            elif "intronic" in annotation._AnnovarParser__Func_refGene:
                export_string += "Tolerated (Intron)"
            elif "Intron" in snp.get_context():
                export_string += "Tolerated (Intron)"
            elif "Intergenic" in snp.get_context():
                export_string += "Tolerated (Intergenic)"
            elif "Coding" in snp.get_context() and "synonymous_variant" in snp.get_consequences():
                export_string += "Tolerated (synonymous_variant)"
            else:
                export_string += "."

    # change digit format to german
    export_string = export_string.replace('.', ',')
    # filter deleterious only
    if args.filter_deleterious and annotation._AnnovarParser__MetaLR_pred == "D":
        target.write(export_string)
        target.write("\n")
    if not args.filter_deleterious:
        target.write(export_string)
        target.write("\n")
    export_cnt += 1
target.close()

# HTML page generation ###
doc = dominate.document(title='SAPA output')
export_cnt = 0
wd = os.path.dirname(os.path.realpath(__file__))

with doc.head:
    # link(rel='stylesheet', href='https://cdn.datatables.net/1.10.13/css/jquery.dataTables.min.css')
    link(rel='stylesheet', href='./resources/dataTables.min.css')
    link(rel='stylesheet', href='./resources/SAPA.css')
    # script(type='text/javascript', src='https://code.jquery.com/jquery-1.12.4.js')
    script(type='text/javascript', src='./resources/jquery.1.12.4.min.js')
    # script(type='text/javascript', src='https://cdn.datatables.net/1.10.13/js/jquery.dataTables.min.js')
    script(type='text/javascript', src='./resources/dataTables.min.js')
    # script(type='text/javascript', src='https://cdn.datatables.net/plug-ins/1.10.13/sorting/natural.js')
    script(type='text/javascript', src='./resources/natural.js')
    script(type='text/javascript', src="./resources/main.js")

with doc.add(div(id='content')):
    h1('SAPA - SNP Annotation Programm for AML')
    info = p('All annotated SNPs are listed below, see ')
    info.add(a("SAPA GIT", href='https://github.com/Twinstar2/SAPA'))
    info.add(" for details.")
    button("Invert colors".title(), id="invert_button")

    with table(id="resultTable", border='1'):
        # line = tr()
        thead = thead()
        tbody = tbody()

    # write header first:
    if args.detail:
        header = str(snps_with_annotation.iterkeys().next().print_header()) + \
                 str(snps_with_annotation.itervalues().next().print_header() + "final prediction\t")
    else:
        header = str(snps_with_annotation.iterkeys().next().print_header()) + \
                 "Gene\tfunction prediction scores [0-1]\tconservation scores[-12.3-6.17]\t" \
                 "ensemble scores[0-60]\tfinal prediction\t"

    split_header = header.split("\t")
    for i in split_header:
        if i != "":
            thead += td(i)
    # thead += tr()
    print "writing header"

    for snp, annotation in ordered_snps_with_annotation.iteritems():
        row = tr()
        # write rows in table
        try:
            altFreq = snp.get_altFreq() * 100
        except ValueError:
            altFreq = snp.get_altFreq()

        if args.detail:
            for value in [snp.get_id(), snp.get_chr(), snp.get_pos(), snp.get_ref(), snp.get_alt(),
                          snp.get_type(), ','.join(snp.get_context()), ','.join(snp.get_consequences()),
                          snp.get_dbSNP(), snp.get_cosmic(), snp.get_clinVar(), snp.get_qual(),
                          altFreq, snp.get_totalDepth(), snp.get_refDepth(),
                          snp.get_altDepth(), snp.get_strandBias(),
                          annotation._AnnovarParser__Chr,
                          annotation._AnnovarParser__Start,
                          annotation._AnnovarParser__End,
                          annotation._AnnovarParser__Ref,
                          annotation._AnnovarParser__Alt,
                          annotation._AnnovarParser__Func_refGene,
                          annotation._AnnovarParser__Gene_refGene,
                          annotation._AnnovarParser__GeneDetail_refGene,
                          annotation._AnnovarParser__ExonicFunc_refGene,
                          annotation._AnnovarParser__AAChange_refGene,
                          annotation._AnnovarParser__cytoBand,
                          annotation._AnnovarParser__esp6500siv2_all,
                          annotation._AnnovarParser__avsnp147,
                          annotation._AnnovarParser__SIFT_score,
                          annotation._AnnovarParser__SIFT_pred,
                          annotation._AnnovarParser__Polyphen2_HDIV_score,
                          annotation._AnnovarParser__Polyphen2_HDIV_pred,
                          annotation._AnnovarParser__Polyphen2_HVAR_score,
                          annotation._AnnovarParser__Polyphen2_HVAR_pred,
                          annotation._AnnovarParser__LRT_score,
                          annotation._AnnovarParser__LRT_pred,
                          annotation._AnnovarParser__MutationTaster_score,
                          annotation._AnnovarParser__MutationTaster_pred,
                          annotation._AnnovarParser__MutationAssessor_score,
                          annotation._AnnovarParser__MutationAssessor_pred,
                          annotation._AnnovarParser__FATHMM_score,
                          annotation._AnnovarParser__FATHMM_pred,
                          annotation._AnnovarParser__PROVEAN_score,
                          annotation._AnnovarParser__PROVEAN_pred,
                          annotation._AnnovarParser__VEST3_score,
                          annotation._AnnovarParser__CADD_raw,
                          annotation._AnnovarParser__CADD_phred,
                          annotation._AnnovarParser__DANN_score,
                          annotation._AnnovarParser__fathmm_MKL_coding_score,
                          annotation._AnnovarParser__fathmm_MKL_coding_pred,
                          annotation._AnnovarParser__MetaSVM_score,
                          annotation._AnnovarParser__MetaSVM_pred,
                          annotation._AnnovarParser__MetaLR_score,
                          annotation._AnnovarParser__MetaLR_pred,
                          annotation._AnnovarParser__integrated_fitCons_score,
                          annotation._AnnovarParser__integrated_confidence_value,
                          annotation._AnnovarParser__GERP_RS,
                          annotation._AnnovarParser__phyloP7way_vertebrate,
                          annotation._AnnovarParser__phyloP20way_mammalian,
                          annotation._AnnovarParser__phastCons7way_vertebrate,
                          annotation._AnnovarParser__phastCons20way_mammalian,
                          annotation._AnnovarParser__SiPhy_29way_logOdds
                          ]:
                if type(value) is tuple and not value[0] == ".":  # (value, min, max, threshold)
                    if len(value) == 3:  # has no threshold
                        range = float(value[2]) - float(value[1])  # max - min
                        percentage = float(value[0]) / range
                        red = 255 * percentage
                        green = 255 - red
                        blue = 20
                        colorstring = str(int(red)) + ", " + str(int(green)) + ", " + str(blue)
                        row += td(value[0], style='background-color: rgb(' + colorstring + ")")
                    elif value == annotation._AnnovarParser__SIFT_score or \
                         value == annotation._AnnovarParser__LRT_score or \
                         value == annotation._AnnovarParser__FATHMM_score or \
                         value == annotation._AnnovarParser__PROVEAN_score:
                        # reverse score !!! ###
                        if float(value[0]) < float(value[3]):  # below threshold
                            # red -> yellow
                            range = float(value[3]) - float(value[1])
                            percentage = float(value[0]) / range
                            red = 255
                            green = 255 * percentage
                            blue = 20
                            colorstring = str(int(red)) + ", " + str(int(green)) + ", " + str(blue)
                            row += td(value[0], style='background-color: rgb(' + colorstring + ")")
                        else:  # over threshold
                            # yellow -> red
                            range = float(value[2]) - float(value[3])
                            percentage = float(value[0]) / range
                            if percentage > 1:
                                percentage = 1
                            red = 255 * (1-percentage)
                            green = 255
                            blue = 20
                            colorstring = str(int(red)) + ", " + str(int(green)) + ", " + str(blue)
                            row += td(value[0], style='background-color: rgb(' + colorstring + ")")
                    else:  # has threshold, normal score ###
                        # print value
                        if float(value[0]) < float(value[3]):  # below threshold
                            # green -> yellow
                            range = float(value[3]) - float(value[1])
                            percentage = float(value[0]) / range
                            red = 255 * percentage
                            green = 255
                            blue = 20
                            colorstring = str(int(red)) + ", " + str(int(green)) + ", " + str(blue)
                            row += td(value[0], style='background-color: rgb(' + colorstring + ")")
                        else:  # over threshold
                            # yellow -> red
                            range = float(value[2]) - float(value[3])  # max - min
                            if range < 0:
                                range = abs(range)
                            percentage = float(value[0]) / range
                            if value == annotation._AnnovarParser__DANN_score:
                                print str(value[2]) +", "+ str(value[3])
                                print value[0]
                                print range
                                print percentage
                            red = 255
                            green = 255 * (1-percentage)
                            blue = 20
                            colorstring = str(int(red)) + ", " + str(int(green)) + ", " + str(blue)
                            row += td(value[0], style='background-color: rgb(' + colorstring + ")")

                elif type(value) is tuple and value[0] == ".":
                    row += td(value[0])
                else:
                    row += td(value)

        else:
            for value in [snp.get_id(), snp.get_chr(), snp.get_pos(),
                          snp.get_ref(), snp.get_alt(), snp.get_type(),
                          ','.join(snp.get_context()), ','.join(snp.get_consequences()),
                          snp.get_dbSNP(), snp.get_cosmic(), snp.get_clinVar(),
                          snp.get_qual(), snp.get_altFreq(),
                          snp.get_totalDepth(), snp.get_refDepth(),
                          snp.get_altDepth(), snp.get_strandBias(),
                          annotation._AnnovarParser__Gene_refGene,
                          annotation._AnnovarParser__MetaLR_score,
                          annotation._AnnovarParser__GERP_RS,
                          annotation._AnnovarParser__CADD_phred
                          ]:

                if type(value) is tuple and not value[0] == ".":
                    range = float(value[2]) - float(value[1])
                    percentage = float(value[0]) / range
                    red = 255 * percentage
                    green = 255 - red
                    blue = 20
                    colorstring = str(int(red)) + ", " + str(int(green)) + ", " + str(blue)
                    #yellow = "r":255, "g":255, "b":20
                    row += td(value[0], style='background-color: rgb(' + colorstring + ")")
                elif type(value) is tuple and value[0] == ".":
                    row += td(value[0])
                else:
                    row += td(value)

        export_string = ""
        # check for final prediction ###
        if annotation._AnnovarParser__MetaLR_pred == "T":
            export_string += "Tolerated"
        elif annotation._AnnovarParser__MetaLR_pred == "D":
            export_string += "Deleterious"
        elif annotation._AnnovarParser__MetaLR_pred == "P":
            export_string += "possibly damaging"
        elif annotation._AnnovarParser__MetaLR_pred == ".":
            # no Meta Socore available
            if annotation._AnnovarParser__fathmm_MKL_coding_pred == "D":
                export_string += "Deleterious (FATHMM)"
            elif not annotation._AnnovarParser__DANN_score[0] == ".":
                if float(annotation._AnnovarParser__DANN_score[0]) > 0.96:
                    export_string += "Deleterious (DANN)"
            else:
                if annotation._AnnovarParser__ExonicFunc_refGene == "synonymous SNV":
                    export_string += "Tolerated (synonymous)"
                elif "intronic" in annotation._AnnovarParser__Func_refGene:
                    export_string += "Tolerated (Intron)"
                elif "Intron" in snp.get_context():
                    export_string += "Tolerated (Intron)"
                elif "Intergenic" in snp.get_context():
                    export_string += "Tolerated (Intergenic)"
                elif "Coding" in snp.get_context() and "synonymous_variant" in snp.get_consequences():
                    export_string += "Tolerated (synonymous_variant)"
                else:
                    export_string += "."

        row += td(export_string)
        tbody += row

        # change digit format to german
        export_string = export_string.replace('.', ',')
        # filter deleterious only
        if args.filter_deleterious and annotation._AnnovarParser__MetaLR_pred == "D":
            # target.write(export_string)
            # target.write("\n")
            print "filtering"
            # if not args.filter_deleterious:
            # target.write(export_string)
            # target.write("\n")
        #    print "not filtering"
        export_cnt += 1

target_html.write(doc.render())
# print(doc.render())
target_html.close()

print "FINISHED"
