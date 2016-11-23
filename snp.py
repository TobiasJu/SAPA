# mutation class for storing the patient data


class SNP(object):
    __id = 0
    __chr = ""
    __pos = 0  # int
    __ref = ""
    __alt = ""  # mutated Base(es)
    __type = ""
    __context = ""
    __consequence = ""
    __dbSNP = ""
    __cosmic = ""
    __clinVar = ""
    __qual = 0  # int
    __altFreq = 0.0  # float
    __totalDepth = 0  # int
    __refDepth = 0  # int
    __altDepth = 0  # int
    __strandBias = 0.0  # float
    __new_end = 0

    # constructor
    def __init__(self, id, chr, pos, ref, alt, type, context, consequence, dbSNP, cosmic, clinVar, qual, altFreq,
                 totalDepth, refDepth, altDepth, strandBias):
        self.__id = id
        self.__chr = chr
        self.__pos = pos
        self.__new_end = pos
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

    # setter
    def set_id(self, id):
        self.__id = id

    def set_chr(self, chr):
        self.__chr = chr

    def set_pos(self, pos):
        self.__pos = pos
        self.__new_end = pos

    def set_ref(self, ref):
        self.__ref = ref

    def set_alt(self, alt):
        self.__alt = alt

    def set_type(self, type):
        self.__type = type

    def set_context(self, context):
        if len(context.split(',')) == 1:
            self.__context = context
        else:
            multiple_context = context.split(',')
            self.__context = multiple_context

    def set_consequence(self, consequence):
        if len(consequence.split(',')) == 1:
            self.__context = consequence
        else:
            multiple_context = consequence.split(',')
            self.__context = multiple_context

    def set_dbSNP(self, dbSNP):
        if len(dbSNP.split(',')) == 1:
            self.__context = dbSNP
        else:
            multiple_context = dbSNP.split(',')
            self.__context = multiple_context

    def set_cosmic(self, cosmic):
        if len(cosmic.split(',')) == 1:
            self.__context = cosmic
        else:
            multiple_context = cosmic.split(',')
            self.__context = multiple_context

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

    def set_new_end(self, new_end):
        self.__new_end = new_end

    # getter
    def get_id(self):
        return int(self.__id)

    def get_chr(self):
        return self.__chr

    def get_pos(self):
        return int(self.__pos)

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
        return int(self.__qual)

    def get_altFreq(self):
        return float(self.__altFreq)

    def get_totalDepth(self):
        return int(self.__totalDepth)

    def get_refDepth(self):
        return int(self.__refDepth)

    def get_altDepth(self):
        return int(self.__altDepth)

    def get_strandBias(self):
        return float(self.__strandBias)

    def get_new_end(self):
        if self.__new_end == 0:
            return self.__pos
        else:
            return self.__new_end

    # def get_conclusion(self):
    #     return self.__conclusion

    # class functions
    def toString(self):
        return "ID: {}, Chr: {}, Pos: {}, Ref: {}, Alt: {}, Type: {}, Context: {}, Consequence: {}, dbSNP: {}, " \
               "COSMIC: {}, ClinVar: {}, Qual: {}, Alt Freq: {}, Total Depth: {}, Ref Depth: {}, Alt Depth: {}, " \
               "Strand Bias: {} ".format(self.__id,
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
                                         self.__strandBias)

    def print_header(self):
        return "ID\tChr\tPos\tRef\tAlt\tType\tContext\tConsequence\tdbSNP\t" \
               "COSMIC\tClinVar\tQual\tAlt Freq\tTotal Depth\tRef Depth\tAlt Depth\t" \
               "Strand Bias\t"

    def export(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\n".format(self.__chr,
                                                                                               self.__pos,
                                                                                               self.__new_end,
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
                                                                                               self.__strandBias)
