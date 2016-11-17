#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import csv
import sys
import os
from os.path import exists
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
parser.add_argument("-i", "--input_file", help="tab separated table with SNP's")
parser.add_argument("-s", "--separator", help='set the input file separator (default: ",")')
parser.add_argument("-t", "--text_delimiter", help='set the input text delimiter (default: ")')
parser.add_argument("-d", "--input_directory", type=str, help="hg19 database directory (default /hg19)")
parser.add_argument("-o", "--output_file", type=str, help="output file name (default output.txt)")
parser.add_argument("-f", "--fast", action="store_true", help="run annotation just with a region based approach, "
                                                              "for faster computing and less download file demand")
parser.add_argument("-fi", "--filter", action="store_true", help="filter SNPs for nonsynonymus and clinically "
                                                                 "significant (>95%%) SNPs")
args = parser.parse_args()

# sanity check ###
if not len(sys.argv) > 1:
    parser.print_help()
    sys.exit(0)

if args.manual:
    print """Mutation Information / SNP Information

This program is designed to add additional information to an Illumina truseq amplicon variants file.

The following Parameters are added to the file:

- COMING SOON

Please Note that during the first run of the program, the required database will be downloaded. This may take some time,
 depending on your internet connection.

Several commonly used databases are integrated: ‘cytoBand’ for the chromosome coordinate of each cytogenetic band,
‘1000g2014oct’ for alternative allele frequency in the 1000 Genomes Project (version October 2014), ‘exac03’
for the variants reported in the Exome Aggregation Consortium (version 0.3)50, ‘ljb26_all’ for various functional
deleteriousness prediction scores from the dbNSFP database (version 2.6)51, ‘clinvar_20140929’ for the variants
reported in the ClinVar database (version 20140929)52 and ‘snp138’ for the dbSNP database (version 138)53.
             """
    parser.print_help()
    sys.exit(0)

if not args.input_file:
    print "ERROR, please enter a input file parameter"
    parser.print_help()
    sys.exit(0)

# if not args.input_directory:
#    print "ERROR, please enter a input directory"
#    parser.print_help()
#    sys.exit(0)

print args

# if args.verbose:
# print "detailed output selected"

# load the input file into a variable
with open(args.input_file) as csvfile:
    if args.separator and args.text_delimiter:
        variant_lines = csv.reader(csvfile, args.separator, args.text_delimiter)
    elif args.separator and not args.text_delimiter or args.text_delimiter and not args.separator:
        print "please enter delimiter AND quote chars"
        sys.exit(0)
    else:
        variant_lines = csv.reader(csvfile, delimiter=',', quotechar='"')

    # create mutation Objects from the given patient Data ###

    # skip header if there
    has_header = csv.Sniffer().has_header(csvfile.read(1024))
    csvfile.seek(0)  # rewind
    incsv = csv.reader(csvfile)
    if has_header:
        next(variant_lines)
        print "skipping header"
    l_count = 0
    snps = []
    for variant_line in variant_lines:
        # print variant_line
        # remove /n form end of line
        # variant_line = snp_entry.strip()

        variant_line[0] = variant_line[0].translate(None, '"')
        variant_line[-1] = variant_line[-1].translate(None, '"')

        if len(variant_line) >= 16:
            context = variant_line[5].split(",")
            consequences = variant_line[6].split(",")
            snps.append(SNP(l_count, variant_line[0], variant_line[1], variant_line[2], variant_line[3],
                            variant_line[4], context, consequences, variant_line[7], variant_line[8], variant_line[9],
                            int(variant_line[10]), variant_line[11], variant_line[12], variant_line[13],
                            variant_line[14], variant_line[15]))
        else:
            print "INVALID DATA (length) in Line {}".format(l_count)
            print variant_line
        l_count += 1
print "created patient SNP objects with " + str(len(snps)) + " unique SNPs\n"
csvfile.close()

# filter all entries in the patient mutation data set ###
if args.filter:
    print "pre filter SNP count: " + str(len(snps))
    i = 0
    for mutation in snps:
        # only clinically relevant quality
        if mutation.get_qual() <= 95:
            del snps[i]

        # if mutation does not change the amino acid, it does not affect the cell (in most cases)
        if "synonymous_variant" in mutation.get_consequences():
            del snps[i]
        i += 1
    print "past filter SNP count: " + str(len(snps))

# write tab delimited file for annovar #
tab_mutations = open('amplicon_variants_tab.csv', 'w')
for mutation in snps:
    # check if type is deletion, correction of the data for annovar
    if "Deletion" in mutation.get_type():
        ##### IST DAS RICHTIG??
        # print mutation.get_ref()
        # print mutation.get_alt()
        newEnd = mutation.get_pos() + (len(mutation.get_ref()) - 2)
        mutation.set_alt("-")
        if len(mutation.get_ref()) > 2:
            mutation.set_ref(0)
        else:
            mutation.set_ref(mutation.get_ref()[1])
        mutation.set_new_end(newEnd)
        # print mutation.get_pos()
        # print newEnd
        # print mutation.get_ref()
        # print mutation.get_alt()
    tab_mutations.write(mutation.export())
tab_mutations.close()
print "created tab delimited file for annovar"

# WORKS JUST UNDER UBUNTU OR THE UBUNTU BASH FOR WINDOWS #
# get annovar databases if needed ###

annotate_variation = "./perl/annotate_variation.pl "
databases = ["-buildver hg19 -downdb -webfrom annovar refGene hg19/",
             "-buildver hg19 -downdb cytoBand hg19/",
             "-buildver hg19 -downdb -webfrom annovar esp6500siv2_all hg19/",
             "-buildver hg19 -downdb -webfrom annovar 1000g2014oct hg19/",
             "-buildver hg19 -downdb -webfrom annovar snp138 hg19/",
             "-buildver hg19 -downdb -webfrom annovar ljb26_all hg19/"
             ]

if args.fast:
    if not os.path.isfile("hg19/hg19_refGene.txt"):
        print "downloading dependencies..."
        p = subprocess.Popen([annotate_variation + databases[0]], shell=True)
        p.communicate()

    if not os.path.isfile("hg19/hg19_snp138.txt"):
        print "downloading dependencies..."
        p = subprocess.Popen([annotate_variation + databases[5]], shell=True)
        p.communicate()
else:
    if not os.path.isfile("hg19/hg19_refGene.txt"):
        print "downloading dependencies..."
        p = subprocess.Popen([annotate_variation + databases[0]], shell=True)
        p.communicate()

    if not os.path.isfile("hg19/hg19_cytoBand.txt"):
        print "downloading dependencies..."
        p = subprocess.Popen([annotate_variation + databases[1]], shell=True)
        p.communicate()

#    if not os.path.isfile("hg19/hg19_genomicSuperDups.txt"):
#        print "downloading dependencies..."
#        p = subprocess.Popen([annotate_variation + databases[2]], shell=True)
#        p.communicate()

    if not os.path.isfile("hg19/hg19_esp6500siv2_all.txt"):
        print "downloading dependencies..."
        p = subprocess.Popen([annotate_variation + databases[2]], shell=True)
        p.communicate()

    #if os.path.isfile("hg19/hg19_1000g2014oct.zip"):
    #    print "ERROR: unzip is not installed in your system. \nPlease manually uncompress the files " \
    #          "(hg19_1000g2014oct.zip) at the hg19 directory, and rename them by adding hg19_ prefix to the file names."
    #    sys.exit(0)

    #if not os.path.isfile("hg19/hg19_ALL.sites.2014_10.txt"):
    #    print "downloading dependencies..."
    #    p = subprocess.Popen([annotate_variation + databases[3]], shell=True)
    #    p.communicate()

    if not os.path.isfile("hg19/hg19_snp138.txt"):
        print "downloading dependencies..."
        p = subprocess.Popen([annotate_variation + databases[4]], shell=True)
        p.communicate()

    if not os.path.isfile("hg19/hg19_ljb26_all.txt"):
        print "downloading dependencies..."
        p = subprocess.Popen([annotate_variation + databases[5]], shell=True)
        p.communicate()

# run Annovar ###
print "running annovar"
annovar_pl = "./perl/table_annovar.pl "
dir_path = os.path.dirname(os.path.realpath(__file__))
# print dir_path

if args.input_directory:
    annovar_database = args.input_directory
else:
    annovar_database = "hg19/"

if args.fast:
    params = "amplicon_variants_tab.csv " + annovar_database + " -buildver hg19 -out myanno -remove -protocol " \
                                                               "refGene,snp138 -operation g,f -nastring ."
else:
    params = "amplicon_variants_tab.csv " + annovar_database + " -buildver hg19 -out myanno -remove -protocol " \
             "refGene,cytoBand,esp6500siv2_all,snp138,ljb26_all -operation g,r,f,f,f -nastring . "
print annovar_pl + params
# 1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,
#./perl/table_annovar.pl amplicon_variants_tab.csv /hg19/ -buildver hg19 -out myanno -remove -protocol refGene,snp138 -operation g,f -nastring .

p = subprocess.Popen([annovar_pl + params], shell=True)
# wait until it's finished
p.communicate()

# parse annovar file ###
annovar = []
l_count = 0
with open('myanno.hg19_multianno.txt', 'r') as annovar_file:
    for row in annovar_file:
        # filter header
        if l_count != 0:
            row = row.strip()
            row = row.split("\t")
            if len(row) == 6:  # fast run
                annovar.append(AnnovarParser(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8],
                                             row[9], row[10], ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".",
                                             ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".",
                                             ".", ".", ".", ".", ".", "."))
            elif len(row) == 38:
                annovar.append(AnnovarParser(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8],
                                             row[9], row[10], row[11], row[12], row[13], row[14], row[15], row[16],
                                             row[17],
                                             row[18], row[19], row[20], row[21], row[22], row[23], row[24], row[25],
                                             row[26], row[27], row[28], row[29], row[30], row[31], row[32], row[33],
                                             row[34], row[35], row[36], row[37]))
            else:
                print "ERROR could not handle this row: "
                print row
        l_count += 1
annovar_file.close()
print "annotated SNPs count: " + str(len(annovar))

# iterate over annovar data and get final scores ###
for data in annovar:
    score = data._AnnovarParser__SIFT_score
    if score != ".":
        rel_score = float(data._AnnovarParser__SIFT_score) / data._AnnovarParser__SIFT_max
        #print rel_score


###### some work to do here ####






### READ EVALUATION PAPER FIRST!!!

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
print "created all human protein objects"

# search after SNP corresponding genes ###
coding_mutations = []
for mutation in snps:
    if "Coding" in mutation.get_consequences():
        print "{} ID: {}, Position: {}".format("unknown mutation", mutation.get_id(), mutation.get_pos())
        for prot in allHumanProteins:
            if prot.get_start() < mutation.get_pos() < prot.get_end() and mutation.get_chr() == \
                    prot.get_chromosome():
                # print "SNP Pos: {}, Ref Gene Start: {},  Ref Gene End: {}".format(mutation.get_pos(),
                #                                                                        prot.get_start(),
                #                                                                        prot.get_end())
                print "found Gene: " + prot.get_gene()
                coding_mutations.append(ProbedMutation(mutation.get_id(), mutation.get_chr(), mutation.get_pos(),
                                                       mutation.get_ref(), mutation.get_alt(), mutation.get_type(),
                                                       mutation.get_context(), mutation.get_consequences(),
                                                       mutation.get_dbSNP(), mutation.get_cosmic(),
                                                       mutation.get_clinVar(), mutation.get_qual(),
                                                       mutation.get_altFreq(), mutation.get_totalDepth(),
                                                       mutation.get_refDepth(), mutation.get_altDepth(),
                                                       mutation.get_strandBias(), str(prot.get_chromosome()),
                                                       prot.get_gene(), prot.get_geneSyn(), prot.get_geneDesc(),
                                                       prot.get_proteinClass(), prot.get_start(), prot.get_end()
                                                       ))
    else:
        coding_mutations.append(ProbedMutation(mutation.get_id(), mutation.get_chr(), mutation.get_pos(),
                                               mutation.get_ref(), mutation.get_alt(), mutation.get_type(),
                                               mutation.get_context(), mutation.get_consequences(),
                                               mutation.get_dbSNP(), mutation.get_cosmic(),
                                               mutation.get_clinVar(), mutation.get_qual(),
                                               mutation.get_altFreq(), mutation.get_totalDepth(),
                                               mutation.get_refDepth(), mutation.get_altDepth(),
                                               mutation.get_strandBias(), ".", ".", ".", ".", ".", ".", "."
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
    target = open("output.txt", 'w')
    print "writing export in: output.txt"
else:
    target = open(args.output_file, 'w')
    print "writing export in: " + args.output_file
export_cnt = 0
for cmuta in coding_mutations:
    # write header first:
    if export_cnt == 0:
        header = str(cmuta.print_header()) + str(annovar[0].print_header() + "\tcombined score\tESM score")
        target.write(header)
        target.write("\n")
    # write rows in table
    export_string = str("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t"
                        "{}\t"
                        "".format(cmuta.get_id(), cmuta.get_chr(), cmuta.get_pos(),
                                  cmuta.get_ref(), cmuta.get_alt(), cmuta.get_type(),
                                  ','.join(cmuta.get_context()), ','.join(cmuta.get_consequences()),
                                  cmuta.get_dbSNP(), cmuta.get_cosmic(), cmuta.get_clinVar(),
                                  cmuta.get_qual(), cmuta.get_altFreq(),
                                  cmuta.get_totalDepth(), cmuta.get_refDepth(),
                                  cmuta.get_altDepth(), cmuta.get_strandBias(), cmuta.get_geneChromosome(),
                                  cmuta.get_gene(), cmuta.get_geneSyn(), cmuta.get_geneDesc(),
                                  cmuta.get_proteinClass(), cmuta.get_geneStart(), cmuta.get_geneEnd()
                                  ))
    export_string += str(annovar[export_cnt].export_tab())
    target.write(export_string)
    final_score = "\t\t"
    target.write(final_score)
    target.write("\n")
    export_cnt += 1
target.close()

#','.join(
print "FIN"
