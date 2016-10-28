from mutation import Mutation


# initialize super class
class ProbedMutation(Mutation):
    __geneChromosome = ""
    __gene = ""
    __geneSyn = ""
    __geneDesc = ""
    __proteinClass = ""
    __geneStart = 0
    __geneEnd = 0
    __conclusion = ""
    __score = 0.0

    def __init__(self, id, chr, pos, ref, alt, type, context, consequence, dbSNP, cosmic, clinVar, qual, altFreq,
                 totalDepth, refDepth, altDepth, strandBias, geneChromosome, gene, geneSyn, geneDesc, proteinClass,
                 geneStart, geneEnd, conclusion, score):
        self.__geneChromosome = geneChromosome,
        self.__gene = gene,
        self.__geneSyn = geneSyn,
        self.__geneDesc = geneDesc,
        self.__proteinClass = proteinClass,
        self.__geneStart = geneStart,
        self.__geneEnd = geneEnd,
        self.__conclusion = conclusion,
        self.__score = score,
        super(ProbedMutation, self).__init__(id, chr, pos, ref, alt, type, context, consequence, dbSNP, cosmic, clinVar,
                                             qual, altFreq, totalDepth, refDepth, altDepth, strandBias)

    # setters for super class parameters
    def set_geneChromosome(self, geneChromosome):
        self.__geneChromosome = geneChromosome

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

    def set_score(self, score):
        self.__score = score

    # getters for the super class parameters

    def get_geneChromosome(self):
        return str(self.__geneChromosome[0])

    def get_gene(self):
        return str(self.__gene[0])

    def get_geneSyn(self):
        return self.__geneSyn

    def get_geneDesc(self):
        return self.__geneDesc

    def get_proteinClass(self):
        return self.__proteinClass

    def get_geneStart(self):
        return int(self.__geneStart[0])

    def get_geneEnd(self):
        return int(self.__geneEnd[0])

    def get_conclusion(self):
        return self.__conclusion

    def get_score(self):
        return self.__score

    # class functions

    def print_header(self):
        return "ID\tChr\tPos\tRef\tAlt\tType\tContext\tConsequence\tdbSNP\tCOSMIC\tClinVar\tQual\tAlt Freq\t" \
               "Total Depth\tRef Depth\tAlt Depth\tStrand Bias\tGene Chromosome\tGene\tGene synonym\tGene description\t" \
               "Protein class\tGene Start\tGene End\tpathogenic Score\tConclusion"

    # warum gehet das nicht?
    def toString(self):
        return "ID: {}\tChr: {}\tPos: {}\tRef: {}\tAlt: {}\tType: {}\tContext: {}\tConsequence: {}\tdbSNP: {}\t" \
               "COSMIC: {}\tClinVar: {}\tQual: {}\tAlt Freq: {}\tTotal Depth: {}\tRef Depth: {}\tAlt Depth: {}\t" \
               "Strand Bias: {}\tGene Chromosome: {}\tGene: {}\tGene synonym: {}\tGene description: {}\t" \
               "Protein class: {}\tGene Start: {}\tGene End: {}\tpathogenic Score: {}\t" \
               "conclusion: {}".format(self.__id,
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
                                       self.__geneChromosome,
                                       self.__gene,
                                       self.__geneSyn,
                                       self.__geneDesc,
                                       self.__proteinClass,
                                       self.__geneStart,
                                       self.__geneEnd,
                                       self.__score,
                                       self.__conclusion)

    # warum gehet das auch nicht?
    def generate_export(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t" \
               "".format(self.__id,
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
                         self.__geneStart,
                         self.__geneEnd,
                         self.__conclusion,
                         self.__score)
