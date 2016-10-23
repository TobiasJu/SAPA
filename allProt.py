# class for all human proteins

class AllProt:
    __gene = ""
    __geneSyn = ""
    __ensembl = ""
    __geneDesc = ""
    __chromosome = ""
    __start = 0
    __end = 0
    __proteinClass = ""
    __REST = ""

    # constructor
    def __init__(self, gene, geneSyn, ensembl, geneDesc, chromosome, start, end, proteinClass, REST):
        self.__gene = gene,
        self.__geneSyn = geneSyn,
        self.__ensembl = ensembl,
        self.__geneDesc = geneDesc,
        self.__chromosome = chromosome,
        self.__start = start,
        self.__end = end,
        self.__proteinClass = proteinClass,
        self.__REST = REST

    # setter
    def set_gene(self, gene):
        self.__gene = gene

    def set_start(self, start):
        self.__start = start

    # more setters needed?

    # getter
    def get_gene(self):
        return self.__gene

    def get_geneSyn(self):
        return self.__geneSyn

    def get_ensembl(self):
        return self.__ensembl

    def get_geneDesc(self):
        return self.__geneDesc

    def get_chromosome(self):
        return self.__chromosome

    def get_start(self):
        # wieso ist das ein TUPEL???

        # print type(self.__start)
        return int(self.__start[0])

    def get_end(self):
        # ist auch ein Tupel
        # print int(self.__end[0])
        return int(self.__end[0])

    def get_proteinClass(self):
        return self.__proteinClass

    # class functions

    def toString(self):
        return "Gene: {}, Gene synonym: {}, Ensembl: {}, Gene description: {}, Chromosome: {}, Start: {}, " \
               "End: {},Protein class: {}".format(self.__gene,
                                                  self.__geneSyn,
                                                  self.__ensembl,
                                                  self.__geneDesc,
                                                  self.__chromosome,
                                                  self.__start,
                                                  self.__end,
                                                  self.__proteinClass)
                                                # self.__REST
