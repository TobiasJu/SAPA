
class GeneDNA:
    __name = ""
    __chromosome = ""
    __start = 0
    __end = 0
    __na_sequence = ""
    __aa_sequence = ""

    # constructor
    def __init__(self, name, chromosome, start, end, na_sequence, aa_sequence):
        self.__name = name
        self.__chromosome = chromosome
        self.__start = start
        self.__end = end
        self.__na_sequence = na_sequence
        self.__aa_sequence = aa_sequence


    def set_name(self, name):
        self.__name = name

    def set_chromosome(self, chrmosome):
        self.__chromosome = chrmosome

    def set_start(self, start):
        self.__start = start

    def set_end(self, end):
        self.__end = end

    def set_na_sequence(self, na_sequence):
        self.__na_sequence = na_sequence

    def set_aa_sequence(self, aa_sequence):
        self.__aa_sequence = aa_sequence

    # getter
    def get_name(self):
        return self.__name

    def get_start(self):
        # print type(self.__start)
        return int(self.__start[0])

    def get_end(self):
        # print int(self.__end[0])
        return int(self.__end[0])

    def get_na_sequence(self):
        return self.__na_sequence

    def get_aa_sequence(self):
        return self.__aa_sequence