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
l_count = 0

# create mutation Objects from the given Data of a patient
mutations = []
for line in f:
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

print "created User Mutation Objects\n"
f.close()
line_count = 0

# create Objects containing all human proteins
allHumanProteins = []
file = open('data/allprots.csv', 'r')
lines = file.readlines()[1:]

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

        geneSyn = split_line[1].split(",")
        position = split_line[5].split("-")
        start = position[0]
        end = position[1]
        geneDesc = split_line[3].split(",")
        #print split_line
        allHumanProteins.append(AllProt(split_line[0], geneSyn, split_line[2], split_line[3], geneDesc, start,
                                        int(end), split_line[4], split_line[5]))

    # else:
        # print len(split_line)
        # print split_line
        # allHumanProteins.append(AllProt(split_line[0], "unknown", "unknown", "unknown", "unknown", "unknown", "unknown",
        #                               "unknown", "unknown"))
# print "{} -> {}".format(start, end)
file.close()
print "created All Protein Objects"
# print(mutations[0].toString)
print ("Mutation count: ", l_count)

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

export_list = []
# search after unknown mutations and find corresponding genes
for mutation in mutations:
    if mutation.get_dbSNP() == "" and mutation.get_cosmic() == "" and mutation.get_clinVar() == "" \
            and "Coding" in mutation.get_context():

        print "{} ID: {}, Position: {}".format("unknown mutation", mutation.get_id(), mutation.get_pos())

        for gene in allHumanProteins:
            #print type(gene.get_start())
            #print type(gene.get_end())
            #print type(mutation.get_pos())
            if gene.get_start() < mutation.get_pos() < gene.get_end():
                print "Mutation Pos: {}, Ref Gene Start: {},  Ref Gene End: {}".format(mutation.get_pos(),
                                                                                       gene.get_start(), gene.get_end())
                print "found Gene: "
                # print gene.get_gene()
                # print gene.get_geneSyn()
                print gene.get_geneDesc()
                # export_list.append(ProbedMutation()

                # FAM83A	BJ-TSA-9, MGC14128	ENSG00000147689	Family with sequence similarity 83, member A	8	123178960-123210079 POS 123195662

    export_list.append(ProbedMutation(mutation[0],mutation[1],mutation[1],mutation[1],mutation[1],mutation[1],
                                      mutation[1],mutation[1],mutation[1],mutation[1],mutation[1],mutation[1],
                                      mutation[1],mutation[1],mutation[1],mutation[1],mutation[1],mutation[1],))

# ensemble API
# ensembl_rest.run(species="human", symbol="BRAF")


# write in export table
target = open("output.csv", 'w')
export_cnt = 0
for mutation in mutations:
    # write Header first:
    if export_cnt == 0:
        target.write(ProbedMutation[0].print_header())
        export_cnt += 1

    target.write(mutation.export())

target.close()
