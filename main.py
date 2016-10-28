import sys
import os
import ensembl_rest
from mutation import Mutation
from allProt import AllProt
from probedMutation import ProbedMutation
import argparse

# argparse for information
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
group.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
group.add_argument("-q", "--quiet", action="store_true", help="prevent output")
parser.add_argument("-o", "--output", action="store_true", help="store the output in a file")
parser.add_argument("-n", "--number", type=int, help="number X")

args = parser.parse_args()

# if args.verbose:
# do something

f = open('data/truseq-amplicon-variants_tobi.csv', 'r')
variant_lines = f.readlines()[1:]
l_count = 0

# create mutation Objects from the given Data of a patient
mutations = []
for line in variant_lines:
    # remove /n form end of line
    line = line.strip()
    split_line = line.split('","')
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
        print "INVALID DATA LENGHT in Line {}".format(l_count)
        print line
        l_count += 1

print "created User Mutation Objects\n"
f.close()
line_count = 0

# create Objects containing all human proteins
allHumanProteins = []
allProtFile = open('data/allprots.csv', 'r')
lines = allProtFile.readlines()[1:]
for line in lines:
    # remove /n form end of line
    line = line.strip()
    split_line = line.split('\t')
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

        # else:
        # print len(split_line)
        # print split_line
        # allHumanProteins.append(AllProt(split_line[0], "unknown", "unknown",
        # "unknown", "unknown", "unknown", "unknown", "unknown", "unknown"))
# print "{} -> {}".format(start, end)
allProtFile.close()
print "created all protein objects"
# print(mutations[0].toString)
print "Mutation count: ", l_count

# filter all entries in the patient mutation data set
i = 0
for mutation in mutations:
    # only clinically relevant quality
    if mutation.get_qual() <= 95:
        del mutations[i]

    # if mutation does not change the amino acid, it does not affect the cell (in nearly all cases)
    if "synonymous_variant" in mutation.get_consequences():
        del mutations[i]
    # more filters to come
    i += 1
print("past filter: ", len(mutations))

coding_mutations = []
# search after unknown mutations and find corresponding genes
for mutation in mutations:
    # print mutation.get_id()
    # print type(mutation.get_context()) = list
    if mutation.get_dbSNP() == "" and mutation.get_cosmic() == "" and mutation.get_clinVar() == "" \
            and "Coding" in mutation.get_context():

        print "{} ID: {}, Position: {}".format("unknown mutation", mutation.get_id(), mutation.get_pos())

        for prot in allHumanProteins:
            # print type(gene.get_gene())
            # print gene.get_gene()
            # print type(gene.get_end())
            # print type(mutation.get_pos())
            if prot.get_start() < mutation.get_pos() < prot.get_end() and mutation.get_chr() == \
                    prot.get_chromosome():
                print "Mutation Pos: {}, Ref Gene Start: {},  Ref Gene End: {}".format(mutation.get_pos(),
                                                                                       prot.get_start(), prot.get_end())
                print "found Gene: "
                print prot.get_gene()
                # print gene.get_geneSyn()
                # rint gene.get_geneDesc()
                # export_list.append(ProbedMutation()

                # FAM83A	BJ-TSA-9, MGC14128	ENSG00000147689	Family with sequence similarity 83, member A	8
                # 123178960-123210079 POS 123195662

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



for cmuta in coding_mutations:
    print cmuta.print_header()
    print cmuta.toString()



# ensemble API
# ensembl_rest.run(species="human", symbol="BRAF")


# write in export table
target = open("output.csv", 'w')
export_cnt = 0
for probed in coding_mutations:
    # write Header first:
    if export_cnt == 0:
        header = str(probed.print_header())
        print type(probed.print_header())
        target.write(header)
        target.write("\n")
        export_cnt += 1

    # print probed.get_id()
    export_string = str("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t"
                        "{}\t{}\t{}\t"
                        "".format(probed.get_id(), probed.get_chr(), probed.get_pos(),
                                  probed.get_ref(), probed.get_alt(), probed.get_type(),
                                  probed.get_context(), probed.get_consequences(),
                                  probed.get_dbSNP(), probed.get_cosmic(), probed.get_clinVar(),
                                  probed.get_qual(), probed.get_altFreq(),
                                  probed.get_totalDepth(), probed.get_refDepth(),
                                  probed.get_altDepth(), probed.get_strandBias(), probed.get_geneChromosome(),
                                  probed.get_gene(), probed.get_geneSyn(), probed.get_geneDesc(),
                                  probed.get_proteinClass(), probed.get_geneStart(), probed.get_geneEnd(),
                                  probed.get_score(), probed.get_conclusion()))

    # print export_string
    # print probed.get_geneStart()
    # print type(probed.generate_export())
    # target.write(probed.generate_export())

    target.write(export_string)
    target.write("\n")

target.close()

print "FIN"
