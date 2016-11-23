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
    __snp138 = ""
    __SIFT_score = 0.0
    __SIFT_pred = ""
    __Polyphen2_HDIV_score = 0.0
    __Polyphen2_HDIV_pred = ""
    __Polyphen2_HVAR_score = 0.0
    __Polyphen2_HVAR_pred = ""
    __LRT_score = 0.0
    __LRT_pred = ""
    __MutationTaster_score = 0.0
    __MutationTaster_pred = ""
    __MutationAssessor_score = 0.0
    __MutationAssessor_pred = ""
    __FATHMM_score = 0.0
    __FATHMM_pred = ""
    __PROVEAN_score = 0.0
    __PROVEAN_pred = 0.0
    __VEST3_score = 0.0
    __CADD_raw = 0.0
    __CADD_phred = 0.0
    __DANN_score = 0.0
    __fathmm_MKL_coding_score = 0.0
    __fathmm_MKL_coding_pred = 0.0
    __MetaSVM_score = 0.0
    __MetaSVM_pred = ""
    __MetaLR_score = 0.0
    __MetaLR_pred = ""
    __integrated_fitCons_score = 0.0
    __integrated_confidence_value = 0.0
    __GERP_RS = 0.0
    __phyloP7way_vertebrate = 0.0
    __phyloP20way_mammalian = 0.0
    __phastCons7way_vertebrate = 0.0
    __phastCons20way_mammalian = 0.0
    __SiPhy_29way_logOdds = 0.0
    # hypothetical max scores
    __SIFT_max = 10.0
    __Polyphen2_HDIV_max = 10.0
    __Polyphen2_HVAR_max = 10.0
    __LRT_max = 10.0
    __MutationTaster_max = 10.0
    __MutationAssessor_max = 10.0
    __FATHMM_max = 10.0
    __RadialSVM_max = 0.0
    __LR_max = 10.0
    __VEST3_max = 10.0

    # constructor
    def __init__(self, Chr, Start, End, Ref, Alt, Func_refGene, Gene_refGene, GeneDetail_refGene, ExonicFunc_refGene,
                 AAChange_refGene, cytoBand, esp6500siv2_all, snp138, SIFT_score, SIFT_pred, Polyphen2_HDIV_score,
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
        self.__snp138 = snp138
        self.__SIFT_score = SIFT_score
        self.__SIFT_pred = SIFT_pred
        self.__Polyphen2_HDIV_score = Polyphen2_HDIV_score
        self.__Polyphen2_HDIV_pred = Polyphen2_HDIV_pred
        self.__Polyphen2_HVAR_score = Polyphen2_HVAR_score
        self.__Polyphen2_HVAR_pred = Polyphen2_HVAR_pred
        self.__LRT_score = LRT_score
        self.__LRT_pred = LRT_pred
        self.__MutationTaster_score = MutationTaster_score
        self.__MutationTaster_pred = MutationTaster_pred
        self.__MutationAssessor_score = MutationAssessor_score
        self.__MutationAssessor_pred = MutationAssessor_pred
        self.__FATHMM_score = FATHMM_score
        self.__FATHMM_pred = FATHMM_pred
        self.__PROVEAN_score = PROVEAN_score
        self.__PROVEAN_pred = PROVEAN_pred
        self.__VEST3_score = VEST3_score
        self.__CADD_raw = CADD_raw
        self.__CADD_phred = CADD_phred
        self.__DANN_score = DANN_score
        self.__fathmm_MKL_coding_score = fathmm_MKL_coding_score
        self.__fathmm_MKL_coding_pred = fathmm_MKL_coding_pred
        self.__MetaSVM_score = MetaSVM_score
        self.__MetaSVM_pred = MetaSVM_pred
        self.__MetaLR_score = MetaLR_score
        self.__MetaLR_pred = MetaLR_pred
        self.__integrated_fitCons_score = integrated_fitCons_score
        self.__integrated_confidence_value = integrated_confidence_value
        self.__GERP_RS = GERP_RS
        self.__phyloP7way_vertebrate = phyloP7way_vertebrate
        self.__phyloP20way_mammalian = phyloP20way_mammalian
        self.__phastCons7way_vertebrate = phastCons7way_vertebrate
        self.__phastCons20way_mammalian = phastCons20way_mammalian
        self.__SiPhy_29way_logOdds = SiPhy_29way_logOdds

        # class functions

    def print_header(self):
        return "Chr\tSNP Start\tSNP End\tRef\tAlt\tFunc_refGene\tGene_refGene\tGeneDetail_refGene\tExonicFunc_refGene" \
               "\tAAChange_refGene\tcytoBand\tesp6500siv2_all\tsnp138\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\t" \
               "Polyphen2_HDIV_pred\t" \
               "Polyphen2_HVAR_score\tPolyphen2_HVAR_pred\tLRT_score\tLRT_pred\tMutationTaster_score\t" \
               "MutationTaster_pred\tMutationAssessor_score\tMutationAssessor_pred\tFATHMM_score\tFATHMM_pred" \
               "\tRadialSVM_score\tRadialSVM_pred\tLR_score\tLR_pred\tVEST3_score\tCADD_raw\tCADD_phred\tGERP_RS\t" \
               "phyloP46way_placental\tphyloP100way_vertebrate\tSiPhy_29way_logOdds\t"

    def export_tab(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t" \
               "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t" \
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
                         self.__snp138,
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
                         self.__RadialSVM_score,
                         self.__RadialSVM_pred,
                         self.__LR_score,
                         self.__LR_pred,
                         self.__VEST3_score,
                         self.__CADD_raw,
                         self.__CADD_phred,
                         self.__GERP_RS,
                         self.__phyloP46way_placental,
                         self.__phyloP100way_vertebrate,
                         self.__SiPhy_29way_logOdds)

#                         self.__g2014oct_all,
#                         self.__g2014oct_afr,
#                         self.__g2014oct_eas,
#                         self.__g2014oct_eur,