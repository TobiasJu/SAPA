#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import csv
import sys
import os
import subprocess
import ensembl_rest
import itertools
from mutation import Mutation
from allProt import AllProt
from probedMutation import ProbedMutation
from dna_to_aa_translator import Translator
from geneDNA import GeneDNA
from annovarParser import AnnovarParser
import argparse

# argparse for information
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
group.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
group.add_argument("-q", "--quiet", action="store_true", help="prevent output")
parser.add_argument("-o", "--output", action="store_true", help="store the output in a file")
parser.add_argument("-n", "--number", type=int, help="just a Test number")
parser.add_argument("-f", "--input_file", help="tab separated table with SNP's")
parser.add_argument("-d", "--input_directory", type=str, help="hg19 directory")

args = parser.parse_args()

# sanity check ###
if not len(sys.argv) > 1:
    parser.print_help()
    sys.exit(0)

if not args.input_file:
    print "ERROR, please enter a input file parameter"
    parser.print_help()
    sys.exit(0)

if not args.input_directory:
    print "ERROR, please enter a input directory"
    parser.print_help()
    sys.exit(0)

print args

# if args.verbose:
# print "detailed output selected"

# load the input file into a variable
with open(args.input_file) as f:
    variant_lines = f.readlines()[1:]

# create mutation Objects from the given patient Data ###
l_count = 0
mutations = []
for dna_line in variant_lines:
    # remove /n form end of line
    dna_line = dna_line.strip()
    split_line = dna_line.split('","')
    split_line[0] = split_line[0].translate(None, '"')
    split_line[-1] = split_line[-1].translate(None, '"')

    if len(split_line) == 16:
        context = split_line[5].split(",")
        consequences = split_line[6].split(",")
        l_count += 1
        mutations.append(Mutation(l_count, split_line[0], split_line[1], split_line[2], split_line[3], split_line[4],
                                  context, consequences, split_line[7], split_line[8], split_line[9],
                                  int(split_line[10]), split_line[11], split_line[12], split_line[13], split_line[14],
                                  split_line[15]))
    else:
        print "INVALID DATA (length) in Line {}".format(l_count)
        print dna_line
        l_count += 1

print "created User Mutation Objects\n"
f.close()

# write tab delimited file for annovar #
tab_mutations = open('amplicon_variants_tab.csv', 'w')
print "creating tab delimited file for annovar"
for mutation in mutations:
    tab_mutations.write(mutation.export())
tab_mutations.close()

# WORKS JUST UNDER UBUNTU OR THE UBUNTU BASH FOR WINDOWS #
# run Annovar ###
print "running Annovar"
dir_path = os.path.dirname(os.path.realpath(__file__))
# print dir_path

params = "amplicon_variants_tab.csv /mnt/c/annovar/humandb/ -buildver hg19 -out myanno -remove -protocol " \
         "refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas," \
         "1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . "  # -csvout

annovar = "./perl/table_annovar.pl "

# p = subprocess.Popen([annovar + params], shell=True)
# wait until it's finished
# p.communicate()

print annovar + params

# parse annovar file ###
annovar = []
l_count = 0
with open('myanno.hg19_multianno.txt', 'r') as annovar_file:
    for row in annovar_file:
        # filter header
        if l_count == 0:
            print "skip header"
        else:
            row.rstrip()
            row = row.split("\t")
            if len(row) == 43:
                annovar.append(AnnovarParser(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8],
                                             row[9], row[10], row[11], row[12], row[13], row[14], row[15], row[16],
                                             row[17],
                                             row[18], row[19], row[20], row[21], row[22], row[23], row[24], row[25],
                                             row[26], row[27], row[28], row[29], row[30], row[31], row[32], row[33],
                                             row[34], row[35], row[36], row[37], row[38], row[39], row[40], row[41],
                                             row[42]))
        l_count += 1
print l_count
annovar_file.close()
print "annotated SNPs"

# iterate over annovar data and get final scores ###
for data in annovar:
    score = data._AnnovarParser__SIFT_score
    if score == ".":
        print "not a number" + score
    else:
        print score
        rel_score = float(data._AnnovarParser__SIFT_score) / data._AnnovarParser__SIFT_max
        print rel_score


### READ EVALUATION PAPER FIRST!!!

sys.exit(0)
# create Objects containing all human proteins ###
allHumanProteins = []
allProtFile = open('data/allprots.csv', 'r')
lines = allProtFile.readlines()[1:]
line_count = 0
for dna_line in lines:
    # remove /n form end of line
    dna_line = dna_line.strip()
    split_line = dna_line.split('\t')
    line_count += 1
    if len(split_line) >= 20:
        i = 0
        for split in split_line:
            split_line[i] = split_line[i].translate(None, "\'")
            i += 1

        prot = split_line[0]
        geneSyn = split_line[1].split(",")
        ensembl = split_line[2]
        position = split_line[5].split("-")
        start = position[0]
        end = position[1]
        geneDesc = split_line[3].split(",")
        chromosome = "chr" + str(split_line[4])
        # print chromosome
        allHumanProteins.append(AllProt(prot, geneSyn, ensembl, geneDesc, chromosome, int(start),
                                        int(end), split_line[6], split_line[7:-1]))

allProtFile.close()
print "created all protein objects"
# print(mutations[0].toString)
print "Mutation count: ", l_count

# filter all entries in the patient mutation data set ###
# i = 0
# for mutation in mutations:
# only clinically relevant quality
#    if mutation.get_qual() <= 95:
#        del mutations[i]

# if mutation does not change the amino acid, it does not affect the cell (in nearly all cases)
#    if "synonymous_variant" in mutation.get_consequences():
#        del mutations[i]
# more filters to come
#    i += 1
# print "past filter: " + str(len(mutations))


# search after unknown mutations and find corresponding genes ###
coding_mutations = []
for mutation in mutations:
    # if mutation.get_dbSNP() == "" and mutation.get_cosmic() == "" and mutation.get_clinVar() == "" \
    # and "Coding" in mutation.get_context():
    if "Coding" in mutation.get_context():
        print "{} ID: {}, Position: {}".format("unknown mutation", mutation.get_id(), mutation.get_pos())
        for prot in allHumanProteins:
            # print type(gene.get_gene())
            # print gene.get_gene()
            # print type(gene.get_end())
            # print type(mutation.get_pos())
            if prot.get_start() < mutation.get_pos() < prot.get_end() and mutation.get_chr() == \
                    prot.get_chromosome():
                # print "Mutation Pos: {}, Ref Gene Start: {},  Ref Gene End: {}".format(mutation.get_pos(),
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
                                                       prot.get_proteinClass(), prot.get_start(), prot.get_end(),
                                                       "non pathogenic", 0.50))

# find DNA sequence for gene in each region and translate it #
mutation_with_sequence = {}
for c_muta in coding_mutations:
    print "Expected gene length: " + str(c_muta.get_geneEnd() - c_muta.get_geneStart())
    openString = args.input_directory + "/chromFa/" + c_muta.get_geneChromosome() + ".fa"
    hg19_chromosome = open(openString, "r")
    with open(openString) as gf:
        chromosome = gf.read()
        chromosome = chromosome.replace(">" + c_muta.get_geneChromosome(), '')
        chromosome = chromosome.replace("\n", '').replace("\r", '').replace("\n\r", '')

    print len(chromosome)
    gene = chromosome[c_muta.get_geneStart():c_muta.get_geneEnd()]
    print len(gene)

    # Translate DNA to AA #
    amino_seq = Translator().translate_dna_sequence(gene)
    # print amino_seq
    g_dna = GeneDNA(c_muta.get_gene(), c_muta.get_geneChromosome(), c_muta.get_geneStart(), c_muta.get_geneEnd(),
                    gene, amino_seq)

    # print "Element name: " + element.get_name()

    mutation_with_sequence[c_muta] = g_dna

for c_muta, g_dna in mutation_with_sequence.iteritems():
    print c_muta.get_gene()
    print c_muta.get_geneChromosome()
    print c_muta.get_pos()
    print g_dna.get_aa_sequence()

# ensemble API for GRCh37 hg19
# ensembl_rest.run(species="human", symbol="DOPEY2")


# write in export table
target = open("output.csv", 'w')
export_cnt = 0
for cmuta in coding_mutations:
    # write Header first:
    if export_cnt == 0:
        header = str(cmuta.print_header())
        # print type(probed.print_header())
        target.write(header)
        target.write("\n")
        export_cnt += 1

    # print probed.get_id()
    export_string = str("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t"
                        "{}\t{}\t{}\t"
                        "".format(cmuta.get_id(), cmuta.get_chr(), cmuta.get_pos(),
                                  cmuta.get_ref(), cmuta.get_alt(), cmuta.get_type(),
                                  cmuta.get_context(), cmuta.get_consequences(),
                                  cmuta.get_dbSNP(), cmuta.get_cosmic(), cmuta.get_clinVar(),
                                  cmuta.get_qual(), cmuta.get_altFreq(),
                                  cmuta.get_totalDepth(), cmuta.get_refDepth(),
                                  cmuta.get_altDepth(), cmuta.get_strandBias(), cmuta.get_geneChromosome(),
                                  cmuta.get_gene(), cmuta.get_geneSyn(), cmuta.get_geneDesc(),
                                  cmuta.get_proteinClass(), cmuta.get_geneStart(), cmuta.get_geneEnd(),
                                  cmuta.get_score(), cmuta.get_conclusion()))

    # print export_string
    # print probed.get_geneStart()
    # print type(probed.generate_export())
    # target.write(probed.generate_export())

    target.write(export_string)
    target.write("\n")

target.close()

print "FIN"
