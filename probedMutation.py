from mutation import Mutation


# initialize super class
class ProbedMutation(Mutation):
    __gene = ""
    __geneSyn = ""
    __geneDesc = ""
    __proteinClass = ""
    __geneStart = 0
    __geneEnd = 0
    __conclusion = ""

    def __init__(self, id, chr, pos, ref, alt, type, context, consequence, dbSNP, cosmic, clinVar, qual, altFreq,
                 totalDepth, refDepth, altDepth, strandBias, gene, geneSyn, geneDesc, proteinClass, geneStart,
                 geneEnd, conclusion):
        self.__gene = gene,
        self.__geneSyn = geneSyn,
        self.__geneDesc = geneDesc,
        self.__proteinClass = proteinClass,
        self.__geneStart = geneStart,
        self.__geneEnd = geneEnd,
        self.__conclusion = conclusion
        super(ProbedMutation, self).__init__(id, chr, pos, ref, alt, type, context, consequence, dbSNP, cosmic, clinVar,
                                             qual, altFreq, totalDepth, refDepth, altDepth, strandBias)

    # setters for super class parameters
    def set_gene(self, gene):
        self.__gene = gene

    def set_geneSyn(self, geneSyn):
        self.__geneSyn = geneSyn

    def set_geneDesc(self, geneDesc):
        self.__geneDesc = geneDesc

    def set_proteinClass(self, proteinClass):
        self.__proteinClass = proteinClass

    def set_geneStart(self, geneStart):
        self.__geneStart = geneStart

    def set_geneEnd(self, geneEnd):
        self.__geneEnd = geneEnd

    def set_conclusion(self, conclusion):
        self.__conclusion = conclusion

    # getters for the super class parameters
    def get_gene(self):
        return self.__gene

    def get_geneSyn(self):
        return self.__geneSyn

    def get_geneDesc(self):
        return self.__geneDesc

    def get_proteinClass(self):
        return self.__proteinClass

    def get_geneStart(self):
        return self.__geneStart

    def get_geneEnd(self):
        return self.__geneEnd

    def get_conclusion(self):
        return self.__conclusion

    def print_header(self):
        return "ID\tChr\tPos\tRef\tAlt\tType\tContext\tConsequence\tdbSNP\tCOSMIC\tClinVar\tQual\tAlt Freq\t" \
               "Total Depth\tRef Depth\tAlt Depth\tStrand Bias\tGene\tGene synonym\tGene description\tProtein class\t" \
               "Gene Start\tGene End\tConclusion"

    # okay warum gehet das nicht?
    def generate_export(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t" \
               "".format(id,
                         chr,
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
                         self.__gene,
                         self.__geneSyn,
                         self.__geneDesc,
                         self.__proteinClass,
                         self.__geneStart,
                         self.__geneEnd,
                         self.__conclusion)
