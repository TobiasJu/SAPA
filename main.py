import sys
import os


class Mutation:
    __id = ""
    __chr = ""
    __pos = 0
    __ref = ""
    __alt = ""  # mutated Base(es)
    __type = ""
    __context = ""
    __consequence = ""
    __dbSNP = ""
    __cosmix = ""
    __clinVar = ""
    __qual = 0.0
    __altFreq = 0.0
    __totalDepth = 0
    __refDepth = 0
    __altDepth = 0
    __strandBias = 0.0

    # constructor
    def __init__(self, chr, id, pos, ref, alt, type, context, consequence, dbSNP, cosmix, clinVar, qual, altFreq,
                 totalDepth, refDepth, altDepth, strandBias):
        self.__id = id
        self.__chr = chr
        self.__pos = pos
        self.__ref = ref
        self.__alt = alt
        self.__type = type
        self.__context = context
        self.__consequence = consequence
        self.__dbSNP = dbSNP
        self.__cosmix = cosmix
        self.__clinVar = clinVar
        self.__qual = qual
        self.__altFreq = altFreq
        self.__totalDepth = totalDepth
        self.__refDepth = refDepth
        self.__altDepth = altDepth
        self.__strandBias = strandBias

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

    def set_cosmix(self, cosmix):
        self.__cosmix = cosmix

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

    


f = open('data/truseq-amplicon-variants_tobi.csv', 'r')
for line in f:
    print line

print "fin\n"
