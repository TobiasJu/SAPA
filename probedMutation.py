from mutation import Mutation

# initialize super class
class ProbedMutation(Mutation):
    __gene = ""
    __geneSyn = ""
    __geneDesc = ""
    __proteinClass = ""
    __conclusion = ""

    def __init__(self, id, chr, pos, ref, alt, type, context, consequence, dbSNP, cosmic, clinVar, qual, altFreq,
                 totalDepth, refDepth, altDepth, strandBias, gene, geneSyn, geneDesc, proteinClass, conclusion):
        super(Mutation, self).__init__(id, chr, pos, ref, alt, type, context, consequence, dbSNP, cosmic, clinVar,
                                       qual, altFreq, totalDepth, refDepth, altDepth, strandBias)
        self.__gene = gene,
        self.__geneSyn = geneSyn,
        self.__geneDesc = geneDesc,
        self.__proteinClass = proteinClass,
        self.__conclusion = conclusion

    # setters for super class parameters
    def set_gene(self, gene):
        self.__gene = gene

    def set_geneSyn(self, geneSyn):
        self.__geneSyn = geneSyn

    def set_geneDesc(self, geneDesc):
        self.__geneDesc = geneDesc

    def set_proteinClass(self, proteinClass):
        self.__proteinClass = proteinClass

    def set_conclusion(self, conclusion):
        self.__conclusion = conclusion

    # getters for the super class parameters
    def get_gene(self):
        return self.__gene

    def get_geneSyn(self):
        return self.__geneSyn

    def get_conclusion(self):
        return self.__conclusion

    def print_header(self):
        return "Chr\tPos\tRef\tAlt\tType\tContext\tConsequence\tdbSNP\tCOSMIC\tClinVar\tQual\tAlt Freq\tTotal Depth\t" \
               "Ref Depth\tAlt Depth\tStrand Bias\tGene\tGene synonym\tGene description\tPosition\tProtein class\t"

    def export(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
            self.__id,
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
            self.__gene,
            self.__geneSyn,
            self.__geneDesc,
            self.__proteinClass,
            self.__conclusion)
