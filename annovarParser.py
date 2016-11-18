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
#    __g2014oct_all = 0.0
#    __g2014oct_afr = 0.0
#    __g2014oct_eas = 0.0
#    __g2014oct_eur = 0.0
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
    __RadialSVM_score = 0.0
    __RadialSVM_pred = ""
    __LR_score = 0.0
    __LR_pred = ""
    __VEST3_score = 0.0
    __CADD_raw = 0.0
    __CADD_phred = 0.0
    __GERP_RS = 0.0
    __phyloP46way_placental = 0.0
    __phyloP100way_vertebrate = 0.0
    __SiPhy_29way_logOdds = 0.0
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
                 AAChange_refGene, cytoBand, esp6500siv2_all, snp138, SIFT_score, SIFT_pred, Polyphen2_HDIV_score, Polyphen2_HDIV_pred,
                 Polyphen2_HVAR_score, Polyphen2_HVAR_pred, LRT_score, LRT_pred, MutationTaster_score,
                 MutationTaster_pred, MutationAssessor_score, MutationAssessor_pred, FATHMM_score, FATHMM_pred,
                 RadialSVM_score, RadialSVM_pred, LR_score, LR_pred, VEST3_score, CADD_raw, CADD_phred,
                 GERP_RS, phyloP46way_placental, phyloP100way_vertebrate, SiPhy_29way_logOdds):
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
        self.__RadialSVM_score = RadialSVM_score
        self.__RadialSVM_pred = RadialSVM_pred
        self.__LR_score = LR_score
        self.__LR_pred = LR_pred
        self.__VEST3_score = VEST3_score
        self.__CADD_raw = CADD_raw
        self.__CADD_phred = CADD_phred
        self.__GERP_RS = GERP_RS
        self.__phyloP46way_placental = phyloP46way_placental
        self.__phyloP100way_vertebrate = phyloP100way_vertebrate
        self.__SiPhy_29way_logOdds = SiPhy_29way_logOdds

        #        self.__g2014oct_all = g2014oct_all
        #        self.__g2014oct_afr = g2014oct_afr
        #        self.__g2014oct_eas = g2014oct_eas
        #        self.__g2014oct_eur = g2014oct_eur

        # g2014oct_all, g2014oct_afr, g2014oct_eas, g2014oct_eur,

    # class functions

    def print_header(self):
        return "Chr\tSNP Start\tSNP End\tRef\tAlt\tFunc_refGene\tGene_refGene\tGeneDetail_refGene\tExonicFunc_refGene" \
               "\tAAChange_refGene\tcytoBand\tesp6500siv2_all\tsnp138\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\t" \
               "Polyphen2_HDIV_pred\t" \
               "Polyphen2_HVAR_score\tPolyphen2_HVAR_pred\tLRT_score\tLRT_pred\tMutationTaster_score\t" \
               "MutationTaster_pred\tMutationAssessor_score\tMutationAssessor_pred\tFATHMM_score\tFATHMM_pred" \
               "\tRadialSVM_score\tRadialSVM_pred\tLR_score\tLR_pred\tVEST3_score\tCADD_raw\tCADD_phred\tGERP_RS\t" \
               "phyloP46way_placental\tphyloP100way_vertebrate\tSiPhy_29way_logOdds\tfunction prediction scores" \
               "\tconservation scores\tensemble scores "

#1000g2014oct_all\t1000g2014oct_afr\t""1000g2014oct_eas" \"\t1000g2014oct_eur\t

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