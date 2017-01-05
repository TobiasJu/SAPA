# annovar class for storing annotated variation data


class AnnovarParser:
    __Chr = ""
    __Start = 0
    __End = 0
    __Ref = ""
    __Alt = ""
    __Func_refGene = ""
    __Gene_refGene = ""
    __GeneDetail_refGene = ""
    __ExonicFunc_refGene = ""
    __AAChange_refGene = ""
    __cytoBand = ""
    __esp6500siv2_all = 0.0
    __avsnp147 = ""
    __SIFT_score = (0.0, 0, 1) # value, min, max
    __SIFT_pred = ""
    __Polyphen2_HDIV_score = (0.0, 0, 1)
    __Polyphen2_HDIV_pred = ""
    __Polyphen2_HVAR_score = (0.0, 0, 1)
    __Polyphen2_HVAR_pred = ""
    __LRT_score = (0.0, 0, 1)
    __LRT_pred = ""
    __MutationTaster_score = (0.0, 0, 1)
    __MutationTaster_pred = ""
    __MutationAssessor_score = (0.0, -5.545, 5.975)
    __MutationAssessor_pred = ""
    __FATHMM_score = (0.0, -16.13, 10.64)
    __FATHMM_pred = ""
    __PROVEAN_score = (0.0, -13, 4)
    __PROVEAN_pred = ""
    __VEST3_score = (0.0, 0, 1)
    __CADD_raw = (0.0, 0, 1)  # ??????????????
    __CADD_phred = (0.0, 0, 60)
    __DANN_score = (0.0, 0, 1)
    __fathmm_MKL_coding_score = (0.0, 0, 1)
    __fathmm_MKL_coding_pred = ""
    __MetaSVM_score = (0.0, -1, 1)
    __MetaSVM_pred = ""
    __MetaLR_score = (0.0, 0, 1.0)
    __MetaLR_pred = ""
    __integrated_fitCons_score = (0.0, 0, 0.8)
    __integrated_confidence_value = 0.0
    __GERP_RS = (0.0, -12.3, 6.17)
    __phyloP7way_vertebrate = (0.0, 0, 1)
    __phyloP20way_mammalian = (0.0, 0, 1)
    __phastCons7way_vertebrate = (0.0, 0, 1)
    __phastCons20way_mammalian = (0.0, 0, 1)
    __SiPhy_29way_logOdds = (0.0, 0, 38)
    # minimal scores
    __SIFT_min = 0
    __Polyphen2_HDIV_min = 0
    __Polyphen2_HVAR_min = 0
    __LRT_min = 0
    __MutationTaster_min = 0
    __MutationAssessor_min = -5.545
    __FATHMM_min = -16.3
    __PROVEAN_min = -13
    __VEST3_min = 0
    __CADD_raw_min = 0
    __CADD_phred_min = 0
    __DANN_min = 0.01
    __fathmm_MKL_coding_min = 0
    __MetaSVM_min = -1
    __MetaLR_min = 0
    __integrated_fitCons_min = 0.0
    __GERP_RS_min = -12.3
    __phyloP7way_vertebrate_min = -1.9
    __phyloP20way_mammalian_min = -3.1
    __phastCons7way_vertebrate_min = 0
    __phastCons20way_mammalian_min = 0
    __SiPhy_29way_logOdds_min = 0
    # observed max scores
    __SIFT_max = 1.0
    __Polyphen2_HDIV_max = 1.0
    __Polyphen2_HVAR_max = 1.0
    __LRT_max = 1.0
    __MutationTaster_max = 1.0
    __MutationAssessor_max = 5.975
    __FATHMM_max = 10.64
    __PROVEAN_max = 4
    __VEST3_max = 1.0
    __CADD_raw_max = 1.0 # ??????
    __CADD_phred_max = 60
    __DANN_max = 0.99
    __fathmm_MKL_coding_max = 1.0
    __MetaSVM_max = 1.0
    __MetaLR_max = 1.0
    __integrated_fitCons_max = 1.0
    __GERP_RS_max = 6.17
    __phyloP7way_vertebrate_max = 1.1
    __phyloP20way_mammalian_max = 1.2
    __phastCons7way_vertebrate_max = 1.0
    __phastCons20way_mammalian_max = 1.0
    __SiPhy_29way_logOdds_max = 38.0


    # constructor
    def __init__(self, Chr, Start, End, Ref, Alt, Func_refGene, Gene_refGene, GeneDetail_refGene, ExonicFunc_refGene,
                 AAChange_refGene, cytoBand, esp6500siv2_all, avsnp147, SIFT_score, SIFT_pred, Polyphen2_HDIV_score,
                 Polyphen2_HDIV_pred, Polyphen2_HVAR_score, Polyphen2_HVAR_pred, LRT_score, LRT_pred,
                 MutationTaster_score, MutationTaster_pred, MutationAssessor_score, MutationAssessor_pred, FATHMM_score,
                 FATHMM_pred, PROVEAN_score, PROVEAN_pred, VEST3_score, CADD_raw, CADD_phred, DANN_score,
                 fathmm_MKL_coding_score, fathmm_MKL_coding_pred, MetaSVM_score, MetaSVM_pred, MetaLR_score,
                 MetaLR_pred, integrated_fitCons_score, integrated_confidence_value, GERP_RS, phyloP7way_vertebrate,
                 phyloP20way_mammalian, phastCons7way_vertebrate, phastCons20way_mammalian, SiPhy_29way_logOdds):
        self.__Chr = Chr
        self.__Start = Start
        self.__End = End
        self.__Ref = Ref
        self.__Alt = Alt
        self.__Func_refGene = Func_refGene
        self.__Gene_refGene = Gene_refGene
        self.__GeneDetail_refGene = GeneDetail_refGene
        self.__ExonicFunc_refGene = ExonicFunc_refGene
        self.__AAChange_refGene = AAChange_refGene
        self.__cytoBand = cytoBand
        self.__esp6500siv2_all = esp6500siv2_all
        self.__avsnp147 = avsnp147
        self.__SIFT_score = (SIFT_score, self.__SIFT_min, self.__SIFT_max, 0.05)  # (value, min, max, Deleterious_threshold)
        self.__SIFT_pred = SIFT_pred
        self.__Polyphen2_HDIV_score = (Polyphen2_HDIV_score, self.__Polyphen2_HDIV_min, self.__Polyphen2_HDIV_max, 0.452)
        self.__Polyphen2_HDIV_pred = Polyphen2_HDIV_pred
        self.__Polyphen2_HVAR_score = (Polyphen2_HVAR_score, self.__Polyphen2_HVAR_min, self.__Polyphen2_HVAR_max, 0.452)
        self.__Polyphen2_HVAR_pred = Polyphen2_HVAR_pred
        self.__LRT_score = (LRT_score, self.__LRT_min, self.__LRT_max, 0.002)
        self.__LRT_pred = LRT_pred
        self.__MutationTaster_score = (MutationTaster_score, self.__MutationTaster_min, self.__MutationTaster_max, 0.5)
        self.__MutationTaster_pred = MutationTaster_pred
        self.__MutationAssessor_score = (MutationAssessor_score, self.__MutationAssessor_min, self.__MutationAssessor_max, 0.65)
        self.__MutationAssessor_pred = MutationAssessor_pred
        self.__FATHMM_score = (FATHMM_score, self.__FATHMM_min, self.__FATHMM_max, -1)
        self.__FATHMM_pred = FATHMM_pred
        self.__PROVEAN_score = (PROVEAN_score, self.__PROVEAN_min, self.__PROVEAN_max, -2.5)
        self.__PROVEAN_pred = PROVEAN_pred
        self.__VEST3_score = (VEST3_score, self.__VEST3_min, self.__VEST3_max)
        self.__CADD_raw = (CADD_raw, self.__CADD_raw_min, self.__CADD_raw_max)
        self.__CADD_phred = (CADD_phred, self.__CADD_phred_min, self.__CADD_phred_max, 15)
        self.__DANN_score = (DANN_score, self.__DANN_min, self.__DANN_max, 0.993)
        self.__fathmm_MKL_coding_score = (fathmm_MKL_coding_score, self.__fathmm_MKL_coding_min, self.__fathmm_MKL_coding_max, 0.63)
        self.__fathmm_MKL_coding_pred = fathmm_MKL_coding_pred
        self.__MetaSVM_score = (MetaSVM_score, self.__MetaSVM_min, self.__MetaSVM_max, 0)
        self.__MetaSVM_pred = MetaSVM_pred
        self.__MetaLR_score = (MetaLR_score, self.__MetaLR_min, self.__MetaLR_max, 0.5)
        self.__MetaLR_pred = MetaLR_pred
        self.__integrated_fitCons_score = (integrated_fitCons_score, self.__integrated_fitCons_min, self.__integrated_fitCons_max)
        self.__integrated_confidence_value = integrated_confidence_value
        self.__GERP_RS = (GERP_RS, self.__GERP_RS_min, self.__GERP_RS_max)
        self.__phyloP7way_vertebrate = (phyloP7way_vertebrate, self.__phyloP7way_vertebrate_min, self.__phyloP7way_vertebrate_max)
        self.__phyloP20way_mammalian = (phyloP20way_mammalian, self.__phyloP20way_mammalian_min, self.__phyloP20way_mammalian_max)
        self.__phastCons7way_vertebrate = (phastCons7way_vertebrate, self.__phastCons7way_vertebrate_min, self.__phastCons7way_vertebrate_max)
        self.__phastCons20way_mammalian = (phastCons20way_mammalian, self.__phastCons20way_mammalian_min, self.__phastCons20way_mammalian_max)
        self.__SiPhy_29way_logOdds = (SiPhy_29way_logOdds, self.__SiPhy_29way_logOdds_min, self.__SiPhy_29way_logOdds_max)

        # class functions

    def print_header(self):
        return "Chr\tStart\tEnd\tRef\tAlt\tFunc_refGene\tGene_refGene\tGeneDetail_refGene\tExonicFunc_refGene\t" \
               "AAChange_refGene\tcytoBand\tesp6500siv2_all\tavsnp147\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\t" \
               "Polyphen2_HDIV_pred\tPolyphen2_HVAR_score\tPolyphen2_HVAR_pred\tLRT_score\tLRT_pred\t" \
               "MutationTaster_score\tMutationTaster_pred\tMutationAssessor_score\tMutationAssessor_pred\tFATHMM_score" \
               "\tFATHMM_pred\tPROVEAN_score\tPROVEAN_pred\tVEST3_score\tCADD_raw\tCADD_phred\tDANN_score\t" \
               "fathmm_MKL_coding_score\tfathmm_MKL_coding_pred\tMetaSVM_score\tMetaSVM_pred\tMetaLR_score\t" \
               "MetaLR_pred\tintegrated_fitCons_score\tintegrated_confidence_value\tGERP_RS\tphyloP7way_vertebrate\t" \
               "phyloP20way_mammalian\tphastCons7way_vertebrate\tphastCons20way_mammalian\tSiPhy_29way_logOdds\t"

    def export_tab(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t" \
               "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t" \
               "".format(self.__Chr,
                         self.__Start,
                         self.__End,
                         self.__Ref,
                         self.__Alt,
                         self.__Func_refGene,
                         self.__Gene_refGene,
                         self.__GeneDetail_refGene,
                         self.__ExonicFunc_refGene,
                         self.__AAChange_refGene,
                         self.__cytoBand,
                         self.__esp6500siv2_all,
                         self.__avsnp147,
                         self.__SIFT_score,
                         self.__SIFT_pred,
                         self.__Polyphen2_HDIV_score,
                         self.__Polyphen2_HDIV_pred,
                         self.__Polyphen2_HVAR_score,
                         self.__Polyphen2_HVAR_pred,
                         self.__LRT_score,
                         self.__LRT_pred,
                         self.__MutationTaster_score,
                         self.__MutationTaster_pred,
                         self.__MutationAssessor_score,
                         self.__MutationAssessor_pred,
                         self.__FATHMM_score,
                         self.__FATHMM_pred,
                         self.__MetaSVM_score,
                         self.__MetaSVM_pred,
                         self.__MetaLR_score,
                         self.__MetaLR_pred,
                         self.__VEST3_score,
                         self.__CADD_raw,
                         self.__CADD_phred,
                         self.__GERP_RS,

                         self.__SiPhy_29way_logOdds)
