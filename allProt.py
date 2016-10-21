
class AllProt:
    __gene = ""
    __geneSyn = ""
    __ensemble = ""
    __geneDesc = ""
    __chromosome = ""
    __position = ""
    __proteinClass = ""
    __REST = ""

    # constructor
    def __init__(self, gene, geneSyn, ensembl, geneDesc, chromosome, position, proteinClass, REST):
        self.__gene = gene,
        self.__geneSyn = geneSyn,
        self.__ensembl = ensembl,
        self.__geneDesc = geneDesc,
        self.__chromosome = chromosome,
        self.__position = position,
        self.__proteinClass = proteinClass,
        self.__REST = REST

    # setter
    def set_gene(self, gene):
        self.__gene = gene

    # more setters needed?

    # getter
    def get_gene(self):
        return self.__gene

    def get_geneSyn(self):
        return self.__geneSyn

    def get_ensembl(self):
        return self.__ensemble

    def get_geneDesc(self):
        return self.__geneDesc

    def get_chromosome(self):
        return self.__chromosome

    def get_position(self):
        return self.__position

    def get_proteinClass(self):
        return self.__proteinClass

    def toString(self):
        return "Gene: {}, Gene synonym: {}, Ensembl: {}, Gene description: {}, Chromosome: {}, Position: {}, " \
               "Protein class: {}".format(self.__gene,
                                          self.__geneSyn,
                                          self.__ensembl,
                                          self.__geneDesc,
                                          self.__chromosome,
                                          self.__position,
                                          self.__proteinClass)
                                         #selt.__REST