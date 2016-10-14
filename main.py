import sys
import os


class Mutation:
    __id = ""
    __chr = ""
    __pos = 0   # int
    __ref = ""
    __alt = ""  # mutated Base(es)
    __type = ""
    __context = ""
    __consequence = ""
    __dbSNP = ""
    __cosmic = ""
    __clinVar = ""
    __qual = 0  # int
    __altFreq = 0.0 # float
    __totalDepth = 0 # int
    __refDepth = 0  # int
    __altDepth = 0  # int
    __strandBias = ""  # 0.0
    __conclusion = ""

    # constructor
    def __init__(self, id, chr, pos, ref, alt, type, context, consequence, dbSNP, cosmic, clinVar, qual, altFreq,
                 totalDepth, refDepth, altDepth, strandBias, conclusion):
        self.__id = id
        self.__chr = chr
        self.__pos = pos
        self.__ref = ref
        self.__alt = alt
        self.__type = type
        self.__context = context
        self.__consequence = consequence
        self.__dbSNP = dbSNP
        self.__cosmic = cosmic
        self.__clinVar = clinVar
        self.__qual = qual
        self.__altFreq = altFreq
        self.__totalDepth = totalDepth
        self.__refDepth = refDepth
        self.__altDepth = altDepth
        self.__strandBias = strandBias
        self.__conclusion = conclusion

    # setter
    def set_id(self, id):
        self.__id = id

    def set_chr(self, chr):
        self.__chr = chr

    def set_pos(self, pos):
        self.__pos = pos

    def set_ref(self, ref):
        self.__ref = ref

    def set_alt(self, alt):
        self.__alt = alt

    def set_type(self, type):
        self.__type = type

    def set_context(self, context):
        self.__context = context

    def set_consequence(self, consequence):
        self.__consequence = consequence

    def set_dbSNP(self, dbSNP):
        self.__dbSNP = dbSNP

    def set_cosmic(self, cosmic):
        self.__cosmic = cosmic

    def set_clinVar(self, clinVar):
        self.__clinVar = clinVar

    def set_qual(selfs, qual):
        selfs.__qual = qual

    def set_altFreq(self, altFreq):
        self.__altFreq = altFreq

    def set_totalDepth(self, totalDepth):
        self.__totalDepth = totalDepth

    def set_refDepth(self, refDepth):
        self.__refDepth = refDepth

    def set_altDepth(self, altDepth):
        self.__altDepth = altDepth

    def set_strandBias(self, strandBias):
        self.__strandBias = strandBias

    # getter
    def get_id(self):
        return self.__id

    def get_chr(self):
        return self.__chr

    def get_pos(self):
        return self.__pos

    def get_ref(self):
        return self.__ref

    def get_alt(self):
        return self.__alt

    def get_type(self):
        return self.__type

    def get_context(self):
        return self.__context

    def get_consequences(self):
        return self.__consequence

    def get_dbSNP(self):
        return self.__dbSNP

    def get_cosmic(self):
        return self.__cosmic

    def get_clinVar(self):
        return self.__clinVar

    def get_qual(self):
        return self.__qual

    def get_altFreq(self):
        return self.__altFreq

    def get_totalDepth(self):
        return self.__totalDepth

    def get_refDepth(self):
        return self.__refDepth

    def get_altDepth(self):
        return self.__altDepth

    def get_strandBias(self):
        return self.__strandBias

    # class functions

    def toString(self):
        return "ID: {}, Chr: {}, Pos: {}, Ref: {}, Alt: {}, Type: {}, Context: {}, Consequence: {}, dbSNP: {}, " \
               "COSMIC: {}, ClinVar: {}, Qual: {}, Alt Freq: {}, Total Depth: {},Ref Depth: {}, Alt Depth: {}, " \
               "Strand Bias: {}, Conclusion: {} ".format(self.__id,
                                         self.__chr,
                                         self.__pos,
                                         self.__ref,
                                         self.__alt,
                                         self.__type,
                                         self.__context,
                                         self.__consequence,
                                         self.__dbSNP,
                                         self.__cosmic,
                                         self.__clinVar,
                                         self.__qual,
                                         self.__altFreq,
                                         self.__totalDepth,
                                         self.__refDepth,
                                         self.__altDepth,
                                         self.__strandBias,
                                         self.__conclusion)


f = open('data/truseq-amplicon-variants_tobi.csv', 'r')
line_count = 0

mutations = []
for line in f:
    # remove /n form end of line
    lien = line.rstrip('/n')
    split_line = line.split('","')
    split_line[0] = split_line[0].translate(None, '"')

    if len(split_line) == 16:
        # wieso geht das nicht?
        # mutation1 = Mutation(split_line)
        line_count += 1
        mutations.append(Mutation(line_count, split_line[0], split_line[1], split_line[2], split_line[3], split_line[4],
                                  split_line[5], split_line[6], split_line[7], split_line[8], split_line[9],
                                  int(split_line[10]), split_line[11], split_line[12], split_line[13], split_line[14],
                                  split_line[15], "unknown"))

print "fin\n"
f.close()

print(mutations[80].toString())
print ("Mutation count: ", line_count)

# filter all entries
i = 0
for mutation in mutations:
    # only clinically relevant quality
    if mutation.get_qual() <= 95:
        del mutations[i]

    # if mutation does not change the amino acid, it does not affect the cell (in nearly all cases)
    if mutation.get_consequences() == "synonymous_variant":
        del mutations[i]
    # more filters to come
    i += 1
print("past filter: ", len(mutations))

for mutation in mutations:
    if mutation.get_clinVar() != "":
        print mutation.get_clinVar()

    if mutation.get_dbSNP() != "":
        print mutation.get_dbSNP()

    if mutation.get_cosmic() != "":
        print mutation.get_cosmic()


# ensebmle API

