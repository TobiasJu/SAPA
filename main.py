import sys
import os
import ensembl_rest
from mutation import Mutation
from allProt import AllProt
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
line_count = 0

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
        line_count += 1
        mutations.append(Mutation(line_count, split_line[0], split_line[1], split_line[2], split_line[3], split_line[4],
                                  context, consequences, split_line[7], split_line[8], split_line[9],
                                  int(split_line[10]), split_line[11], split_line[12], split_line[13], split_line[14],
                                  split_line[15], "unknown"))

print "fin\n"
f.close()
line_count = 0

# create Objects containing all human proteins
allHumanProteins = []
file = open('data/allprots.csv', 'r')

for line in file:
    # remove /n form end of line
    line = line.strip()
    split_line = line.split('\t')
    line_count += 1
    if len(split_line) >= 20:
        allHumanProteins.append(AllProt(split_line[0], split_line[1], split_line[2], split_line[3], split_line[4],
                                        split_line[5], split_line[6], split_line[7]))
    else:
        #print len(split_line)
        #print split_line
        allHumanProteins.append(AllProt(split_line[0], "unknown", "unknown", "unknown", "unknown", "unknown", "unknown",
                                        "unknown"))

file.close()

print allHumanProteins[1].get_geneSyn()

#print(mutations[0].toString)
print ("Mutation count: ", line_count)

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

# search after unknown mutations
for mutation in mutations:
    if mutation.get_dbSNP() == "" and mutation.get_cosmic() == "" and mutation.get_clinVar() == "" \
            and "Coding" in mutation.get_context():

        print "unknown mutation"
        print mutation.get_id()

        # if mutation.get_dbSNP() != "":
        # print mutation.get_dbSNP()

        # if mutation.get_cosmic() != "":
        # print mutation.get_cosmic()


# ensemble API
# ensembl_rest.run(species="human", symbol="BRAF")

# get whole gene, chromosome 
# Suchvorgang...
