SAPA - SNP annotation programm for AML
===================

**75%** of a **Acute myeloid leukemia** (AML) patients SNPs are unique to them and their impact is still unresolved. 

SAPA is designed to add additional information to an **Illumina truseq amplicon** variants csv file. Such as Scores for nonsynonymous scoring matrices, divided into **3** different **combined Scores**, function prediction scores, conservation scores and one ensemble score. So that the viewer gets more Information whether a mutation is deleterious or tolerated. Due to performance and data reasons, SAPA uses annovar<sup>1</sup> and their precomputed SNP database, to quickly gather the deleterious prediction methods scores. 
  

First run 
-------------
<i class="icon-refresh"></i> To get startet simply type in:

`./main.py -i data/truseq_example_data.csv`

this will run the program with the example dataset. Your output will be saved in the `output.txt ` file.  
  

> **Please Note:** During the first run of the program, the required **database will be downloaded**. This may take some time, depending on your internet connection. You may want to use the `-fast` argument to minimize the download size.


Parameter
-------------
![Markdown spport in GIT...](/images/parameters.png)

`-h, --help`
: show the help message and exit

`-q, --quiet`
: prevent output in command line, run the program and don't bother

`-m, --manual`
: display the manual for this program

`-i, --input_file`
: truseq amplicon table containing SNPs 
	> e.g. `-i data/truseq_example_data.csv`

`-d, --detail`
: write detailed output file, with all the single scores of the deleteriousness prediction methods for nonsynonymous SNVs and amino acid substitution (if any) 

`-s, --separator`
: set the input file separator (default: `","`)
	> e.g. `-s ";"`
	this will set the column separator to a semicolon

`-t, --text_delimiter`
: set the input text delimiter (default: `"`)
	> e.g. `-t "'"`
	this will set the text delimiter to an apostrophe, if a column contains multiple entries 

`-b, --buildver` 
: buildversion e.g hg19 (GRCh37) or hg38 (GRCh38) database (default hg19)
	> e.g. `-b hg38`
	set the reference genome to hg38 

`-o, --output_file`
: output file name (default output.txt)
	> e.g. `-o annotated_snps_detailed.txt -d`
	running the detailed operation mode, and saving it in the file annotated_snps_detailed.txt 

`-f, --fast`
: run annotation just with a region based approach, for faster computing and less download file demand

`-fi, --filter`
: filter all entries and only use nonsynonymus and clinically significant (>95%) SNPs in the output file



Features
-------------
#### <i class="icon-pencil"></i>Added data to a SNP csv file

No worries your input file will not be overwritten, a new file with the data from the input file will be created with the following scores added.:

- function prediction scores
- conservation scores
- ensemble score
- final prediction

#### <i class="icon-pencil"></i> Scoring explained (dbnsfp30a)

- **SIFT**
	- Score	ranges from 1 to 0
	- Lower scores are more deleterious
	- The amino acid substitution is predicted damaging if the score is <= 0.05, and tolerated if the score is > 0.05.
	- http://sift.jcvi.org/www/SIFT_help.html
	
- **Polyphen2**
	- Scores ranges from 0 to 1
	- D: Probably damaging (>=0.957), P: possibly damaging (0.453<=pp2<=0.956); B: benign (pp2<=0.452)
	- http://genetics.bwh.harvard.edu/pph2/dokuwiki/overview

- **LRT** 
	- Scores ranges from 1 to 0
	- Deleterious threshold < 0.002 (observed)
	- Lower scores are more deleterious
	- D: Deleterious; N: Neutral; U: Unknown
	- https://www.researchgate.net/profile/Justin_Fay/publication/26671245_Identification_of_deleterious_mutations_within_three_human_genomes/links/00b49529eabbc6aa2b000000.pdf

- **MutationTaster**
	- Scores ranges from 0 to 1
	- Deleterious threshold > 0.5 (https://wikis.utexas.edu/display/bioiteam/Annovar+Annotations)
	- 	A = "disease_causing_automatic"; D =“disease_causing”; N = “polymorphism”; P = “polymorphism_automatic”
	- http://www.mutationtaster.org/


- **MutationAssessor**
	- Scores ranges from -5.545 to 5.975
	- Deleterious threshold > 0.65
	- H: high; M: medium; L: low; N: neutral. H/M means functional and L/N means non-functional
	- http://mutationassessor.org/r3/


- **FATHMM**
	- Scores ranges from -16.13 to 10.64
	- Deleterious threshold < -1 (observed)
	- D: Deleterious; T: Tolerated
	- http://fathmm.biocompute.org.uk/


- **PROVEAN**
	- Scores ranges from -13 to 4
	- Deleterious PROVEAN < -2,5 ; Tolerated PROVEAN > -2,5 
	- http://provean.jcvi.org/about.php


- **VEST3**
	- Scores ranges from 0 to 1
	- The larger the score the more likely the mutation may cause functional change.
	- http://karchinlab.org/apps/appVest.html


- **CADD_raw**
	- Score ranges from -2,375 to 7,6 (observed)
	- The larger the score the more likely the SNP has damaging effect.
	- http://cadd.gs.washington.edu/

- **CADD_phred**
	- Scores range from 0 to 60
	- Deleterious threshold > 15
	- This is [phred-like rank score](https://en.wikipedia.org/wiki/Phred_quality_score) based on whole genome CADD raw scores.
	- The larger the score the more likely the SNP has damaging effect.

- **DANN**
	- Scores ranges from 0.01 to 0.999
	- Deleterious threshold > 0.993 (observed)
	- The higher the score the more likely the mutation may cause functional change.
	- DANN uses the same feature set and training data as CADD to train a deep neural network (DNN)
	- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4341060/

- **FATHMM_MKL_coding_score**
	- Scores ranges from 0 - 1
	- Deleterious threshold > 0.63 (observed)
	- The higher the score, the greater the confidence of a functional mutation.
	- http://fathmm.biocompute.org.uk/fathmmMKL.htm

- **MetaSVM**
	- Scores ranges from -1 to +1
	- Deleterious threshold > 0 (observed)
	- D: Deleterious; T: Tolerated
	- Uses a radial SVM model to train prediction model, using all available scoring algorithm normalized scores
	- http://annovar.openbioinformatics.org/en/latest/user-guide/filter/#-metasvm-annotation

- **MetaLR**
	- Scores ranges from 0 to 1
	- Deleterious threshold > 0.5
	- D: Deleterious; T: Tolerated
	- very similar to MetaSVM, with similar performance
	- http://annovar.openbioinformatics.org/en/latest/user-guide/filter/#-metalr-annotation

- **integrated fitCons**
	- Scores ranges from 0.0 to 1.0
	- higher scores indicating more potential
	    - the fraction of genomic positions evincing a particular pattern (or "fingerprint") of functional assay results, that are under selective pressure
	- http://compgen.cshl.edu/fitCons/

- **GERP RS**
	- RS = (rejected substitutions) 
	- Scores ranges from -12.3 to 6.17
	- The larger the score, the more conserved the site.
	- http://mendel.stanford.edu/SidowLab/downloads/gerp/

- **phyloP7way_vertebrate**
    - Scores ranges from -1.844 to 1.062 (observed)
	- The larger the score, the more conserved.
		- higher scores are more deleterious
	- Phylogenetic p-values for 7 vertebrate species
	- https://brb.nci.nih.gov/seqtools/colexpanno.html#dbnsfp

- **phyloP20way_mammalian**
    - Scores ranges from -3.095 to 1.186 (observed)
	- The larger the score, the more conserved.
		- higher scores are more deleterious
	- Phylogenetic p-values for 20 mammalian species
	- https://brb.nci.nih.gov/seqtools/colexpanno.html#dbnsfp

- **phastCons7way**
	- Scores ranges from 0 to 1
	- The larger the score, the more conserved the site
		- higher scores are more deleterious
	- http://compgen.cshl.edu/phast/phastCons-HOWTO.html

- **phastCons20way**
	- Scores ranges from 0 to 1
	- The larger the score, the more conserved the site
		- higher scores are more deleterious
	- http://compgen.cshl.edu/phast/phastCons-HOWTO.html

- **SiPhy**
	- Scores ranges from 0.0003 to 37.9718
	- The larger the score, the more conserved the site
		- higher scores are more deleterious
	- http://garberlab.umassmed.edu/software.html
                                                                             

### <i class="icon-hdd"></i> Databases
All required databases will be downloaded in the hg19/ <i class="icon-folder-open"></i> directory. 
If some error during downloading occurs, you can read furhter details in the `annovar_downdb.log` in the hg19/ folder.
Several commonly used databases are integrated: 

- <b>cytoBand</b> 
	- for the chromosome coordinate of each cytogenetic band
- <b>exac03</b> 
	- for the variants reported in the Exome Aggregation Consortium (version 0.3)
- <b>dbnsfp30a</b> 
	- for various functional deleteriousness prediction scores from the dbNSFP database (version 2.6)
	
- <b>clinvar_20140929</b>
	- for the variants reported in the ClinVar database (version 20140929).

#### For functional prediction of variants in whole-genome data:

- **GERP++**
	- functional prediction scores for 9 billion mutations based on selective constraints across human genome. You can optionally use gerp++gt2 instead since it includes only RS score greater than 2, which provides high sensitivity while still strongly enriching for truly constrained sites
- **CADD**
	- Combined Annotation Dependent Depletion score for 9 billion mutations. It is basically constructed by a support vector machine trained to differentiate 14.7 million high-frequency human-derived alleles from 14.7 million simulated variants, using ~70 different features.

- **DANN**
	- functional prediction score generated by deep learning, using the identical set of training data as cadd but with much improved performance than CADD.
- **FATHMM**
	- a hidden markov model to predict the functional importance of both coding and non-coding variants (that is, two separate scores are provided) on 9 billion mutations.
- **EIGEN**
	- a spectral approach integrating functional genomic annotations for coding and noncoding variants on 9 billion mutations, without labelled training data (that is, unsupervised approach)
- **GWAVA**
	- genome-wide annotation of variants that supports prioritization of noncoding variants by integrating various genomic and epigenomic annotations on 9 billion mutations.

#### For functional prediction of variants in whole-exome data:

- **dbnsfp30a**
	- this dataset already includes SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, DANN, fitCons, PhyloP and SiPhy scores, but <b>ONLY</b> on coding variants
	- for more information check: http://varianttools.sourceforge.net/Annotation/DbNSFP


#### <i class="icon-file"></i> Input File
As input file is a Illumina truseq amplicon variants file needed, which contains the patients SNPs in a character separated format. A example file is included in the /data directory and listed below.

| Pos       | Ref | Alt | Type     | Context                                    | Consequence                                                                                         | dbSNP                | COSMIC      | ClinVar        | Qual | Alt Freq | Total Depth | Ref Depth | Alt Depth | Strand Bias |
|-----------|-----|-----|----------|--------------------------------------------|-----------------------------------------------------------------------------------------------------|----------------------|-------------|----------------|------|----------|-------------|-----------|-----------|-------------|
| 25463566  | C   | A   | SNV      | Coding                                     | missense_variant                                                                                    |                      |             |                | 100  | 0.334    | 6568        | 4358      | 2194      | -100.0      |
| 128204951 | C   | T   | SNV      | Coding,Intergenic                          | missense_variant,upstream_gene_variant                                                              | rs2335052            | COSM445531  |                | 100  | 0.549    | 1312        | 587       | 720       | -100.0      |
| 128205860 | G   | C   | SNV      | Coding,Intergenic                          | synonymous_variant,upstream_gene_variant                                                            | rs1573858            |             |                | 100  | 0.997    | 1925        | 0         | 1920      | -100.0      |
| 55141055  | A   | G   | SNV      | Coding                                     | synonymous_variant                                                                                  | rs1873778            | COSM1430082 |                | 100  | 0.997    | 12713       | 35        | 12670     | -100.0      |
| 106196829 | T   | G   | SNV      | Coding                                     | missense_variant                                                                                    | rs34402524           | COSM87176   |                | 100  | 0.482    | 38535       | 19801     | 18585     | -100.0      |
| 106196937 | AT  | A   | Deletion | Coding                                     | frameshift_variant,feature_truncation                                                               |                      |             |                | 100  | 0.385    | 31380       | 19297     | 12083     | -100.0      |
| 101917521 | G   | A   | SNV      |                                            |                                                                                                     | rs803064             |             |                | 100  | 0.558    | 9023        | 3973      | 5033      | -100.0      |
| 101921289 | A   | G   | SNV      |                                            |                                                                                                     | rs2230103,rs79334660 |             |                | 100  | 0.37     | 2903        | 1830      | 1073      | -100.0      |
| 21971127  | A   | G   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,synonymous_variant                    |                      |             |                | 100  | 0.172    | 157         | 120       | 27        | -100.0      |
| 21971127  | A   | T   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,synonymous_variant                    |                      |             |                | 44   | 0.057    | 157         | 120       | 9         | -100.0      |
| 21971129  | T   | G   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,missense_variant                      |                      |             |                | 44   | 0.057    | 158         | 143       | 9         | -100.0      |
| 21971130  | G   | A   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,synonymous_variant                    |                      |             |                | 68   | 0.074    | 163         | 144       | 12        | -100.0      |
| 21971131  | G   | A   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,missense_variant                      |                      | COSM13766   |                | 44   | 0.057    | 157         | 144       | 9         | -100.0      |
| 7579472   | G   | C   | SNV      | Coding                                     | missense_variant                                                                                    | rs1042522            | COSM45985   | non-pathogenic | 100  | 0.591    | 3567        | 1452      | 2108      | -100.0      |
| 74733099  | G   | A   | SNV      | Intergenic,Coding,Intron,5P_UTR | downstream_gene_variant,upstream_gene_variant,synonymous_variant,intron_variant,5_prime_UTR_variant | rs237057             |             |                | 100  | 0.995    | 7558        | 34        | 7518      | -100.0      |
| 31022959  | T   | C   | SNV      | Coding                                     | missense_variant                                                                                    | rs6058694            |             |                | 100  | 0.997    | 17877       | 39        | 17831     | -100.0      |
| 36252877  | C   | T   | SNV      | Coding                                     | missense_variant                                                                                    |                      | COSM96546   |                | 100  | 0.465    | 8562        | 4573      | 3978      | -100.0      |
| 39933339  | A   | G   | SNV      | Coding                                     | synonymous_variant                                                                                  | rs5917933            |             |                | 100  | 0.998    | 14616       | 18        | 14592     | -100.0      |
| 123195650 | A   | T   | SNV      | Coding                                     | missense_variant                                                                                    |                      |             |                | 47   | 0.073    | 109         | 100       | 8         | -100.0      |

### <i class="icon-file"></i> Output File (default)
The generated outputfile `output.txt`, which contains the patients annotated SNPs in a character separated format, will be in the current working directory. Using the included example dataset in the /data directory, the following output will be generated.

| ID | Chr   | Pos       | Ref | Alt | Type     | Context                                    | Consequence                                                                                         | dbSNP                | COSMIC      | ClinVar        | Qual | Alt Freq [%] | Total Depth | Ref Depth | Alt Depth | Strand Bias | Gene   | function prediction scores [0-1] | conservation scores[-12.3-6.17] | ensemble scores[0-60] | final prediction       |
|----|-------|-----------|-----|-----|----------|--------------------------------------------|-----------------------------------------------------------------------------------------------------|----------------------|-------------|----------------|------|--------------|-------------|-----------|-----------|-------------|--------|----------------------------------|---------------------------------|-----------------------|------------------------|
| 0  | chr2  | 25463566  | C   | A   | SNV      | Coding                                     | missense_variant                                                                                    |                      |             |                | 100  | 0,334        | 6568        | 4358      | 2194      | -100        | DNMT3A | 0,863                            | 5,76                            | 34                    | Deleterious            |
| 1  | chr3  | 128204951 | C   | T   | SNV      | Coding,Intergenic                          | missense_variant,upstream_gene_variant                                                              | rs2335052            | COSM445531  |                | 100  | 0,549        | 1312        | 587       | 720       | -100        | GATA2  | 0                                | 1,83                            | 13,22                 | Tolerated              |
| 2  | chr3  | 128205860 | G   | C   | SNV      | Coding,Intergenic                          | synonymous_variant,upstream_gene_variant                                                            | rs1573858            |             |                | 100  | 0,997        | 1925        | 0         | 1920      | -100        | GATA2  | 0                                | 0                               | 0                     | Tolerated (synonymous) |
| 3  | chr4  | 55141055  | A   | G   | SNV      | Coding                                     | synonymous_variant                                                                                  | rs1873778            | COSM1430082 |                | 100  | 0,997        | 12713       | 35        | 12670     | -100        | PDGFRA | 0                                | 0                               | 0                     | Tolerated (synonymous) |
| 4  | chr4  | 106196829 | T   | G   | SNV      | Coding                                     | missense_variant                                                                                    | rs34402524           | COSM87176   |                | 100  | 0,482        | 38535       | 19801     | 18585     | -100        | TET2   | 0                                | 5,16                            | 21,5                  | Tolerated              |
| 5  | chr4  | 106196937 | T   | -0  | Deletion | Coding                                     | frameshift_variant,feature_truncation                                                               |                      |             |                | 100  | 0,385        | 31380       | 19297     | 12083     | -100        | TET2   | 0                                | 0                               | 0                     | 0                      |
| 6  | chr7  | 101917521 | G   | A   | SNV      |                                            |                                                                                                     | rs803064             |             |                | 100  | 0,558        | 9023        | 3973      | 5033      | -100        | CUX1   | 0                                | -0,106                          | 5,82                  | Tolerated              |
| 7  | chr7  | 101921289 | A   | G   | SNV      |                                            |                                                                                                     | rs2230103,rs79334660 |             |                | 100  | 0,37         | 2903        | 1830      | 1073      | -100        | CUX1   | 0,004                            | 3,4                             | 7,655                 | Tolerated              |
| 8  | chr9  | 21971127  | A   | G   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,synonymous_variant                    |                      |             |                | 100  | 0,172        | 157         | 120       | 27        | -100        | CDKN2A | 0,192                            | 0,609                           | 1,37                  | Tolerated              |
| 9  | chr9  | 21971127  | A   | T   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,synonymous_variant                    |                      |             |                | 44   | 0,057        | 157         | 120       | 9         | -100        | CDKN2A | 0,178                            | 0,609                           | 0,119                 | Tolerated              |
| 10 | chr9  | 21971129  | T   | G   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,missense_variant                      |                      |             |                | 44   | 0,057        | 158         | 143       | 9         | -100        | CDKN2A | 0,379                            | 5,79                            | 23,7                  | Tolerated              |
| 11 | chr9  | 21971130  | G   | A   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,synonymous_variant                    |                      |             |                | 68   | 0,074        | 163         | 144       | 12        | -100        | CDKN2A | 0,224                            | 3,74                            | 10,53                 | Tolerated              |
| 13 | chr17 | 7579472   | G   | C   | SNV      | Coding                                     | missense_variant                                                                                    | rs1042522            | COSM45985   | non-pathogenic | 100  | 0,591        | 3567        | 1452      | 2108      | -100        | TP53   | 0                                | 1,87                            | 0,355                 | Tolerated              |
| 12 | chr9  | 21971131  | G   | A   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,missense_variant                      |                      | COSM13766   |                | 44   | 0,057        | 157         | 144       | 9         | -100        | CDKN2A | 0,109                            | -0,271                          | 10,15                 | Tolerated              |
| 14 | chr17 | 74733099  | G   | A   | SNV      | Intergenic,Intergenic,Coding,Intron,5P_UTR | downstream_gene_variant,upstream_gene_variant,synonymous_variant,intron_variant,5_prime_UTR_variant | rs237057             |             |                | 100  | 0,995        | 7558        | 34        | 7518      | -100        | SRSF2  | 0                                | 0                               | 0                     | Tolerated (synonymous) |
| 15 | chr20 | 31022959  | T   | C   | SNV      | Coding                                     | missense_variant                                                                                    | rs6058694            |             |                | 100  | 0,997        | 17877       | 39        | 17831     | -100        | ASXL1  | 0                                | 0                               | 0                     | 0                      |
| 16 | chr21 | 36252877  | C   | T   | SNV      | Coding                                     | missense_variant                                                                                    |                      | COSM96546   |                | 100  | 0,465        | 8562        | 4573      | 3978      | -100        | RUNX1  | 0,991                            | 5,31                            | 33                    | Deleterious            |
| 17 | chrX  | 39933339  | A   | G   | SNV      | Coding                                     | synonymous_variant                                                                                  | rs5917933            |             |                | 100  | 0,998        | 14616       | 18        | 14592     | -100        | BCOR   | 0                                | 0                               | 0                     | Tolerated (synonymous) |
| 18 | chrX  | 123195650 | A   | T   | SNV      | Coding                                     | missense_variant                                                                                    |                      |             |                | 47   | 0,073        | 109         | 100       | 8         | -100        | STAG2  | 0,303                            | 4,77                            | 28,6                  | Tolerated              |

### <i class="icon-file"></i> Output File (detailed)
| ID | Chr   | Pos       | Ref | Alt | Type     | Context                                    | Consequence                                                                                         | dbSNP                | COSMIC      | ClinVar        | Qual | Alt Freq [%] | Total Depth | Ref Depth | Alt Depth | Strand Bias | Chr   | Start     | End       | Ref | Alt | Func_refGene | Gene_refGene | GeneDetail_refGene | ExonicFunc_refGene  | AAChange_refGene                                                                                                                                                                                                                                                                                                                                       | cytoBand | esp6500siv2_all | avsnp147    | SIFT_score | SIFT_pred | Polyphen2_HDIV_score | Polyphen2_HDIV_pred | Polyphen2_HVAR_score | Polyphen2_HVAR_pred | LRT_score | LRT_pred | MutationTaster_score | MutationTaster_pred | MutationAssessor_score | MutationAssessor_pred | FATHMM_score | FATHMM_pred | PROVEAN_score | PROVEAN_pred | VEST3_score | CADD_raw | CADD_phred | DANN_score | fathmm_MKL_coding_score | fathmm_MKL_coding_pred | MetaSVM_score | MetaSVM_pred | MetaLR_score | MetaLR_pred | integrated_fitCons_score | integrated_confidence_value | GERP_RS | phyloP7way_vertebrate | phyloP20way_mammalian | phastCons7way_vertebrate | phastCons20way_mammalian | SiPhy_29way_logOdds | final prediction       |
|----|-------|-----------|-----|-----|----------|--------------------------------------------|-----------------------------------------------------------------------------------------------------|----------------------|-------------|----------------|------|--------------|-------------|-----------|-----------|-------------|-------|-----------|-----------|-----|-----|--------------|--------------|--------------------|---------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------|-----------------|-------------|------------|-----------|----------------------|---------------------|----------------------|---------------------|-----------|----------|----------------------|---------------------|------------------------|-----------------------|--------------|-------------|---------------|--------------|-------------|----------|------------|------------|-------------------------|------------------------|---------------|--------------|--------------|-------------|--------------------------|-----------------------------|---------|-----------------------|-----------------------|--------------------------|--------------------------|---------------------|------------------------|
| 0  | chr2  | 25463566  | C   | A   | SNV      | Coding                                     | missense_variant                                                                                    |                      |             |                | 100  | 33,4         | 6568        | 4358      | 2194      | -100        | chr2  | 25463566  | 25463566  | C   | A   | exonic       | DNMT3A       | 0                  | nonsynonymous SNV   | DNMT3A:NM_153759:exon14:c,G1549T:p,G517W,DNMT3A:NM_022552:exon18:c,G2116T:p,G706W,DNMT3A:NM_175629:exon18:c,G2116T:p,G706W                                                                                                                                                                                                                             | 2p23,3   | 0               | rs749365376 | 0          | D         | 1                    | D                   | 1                    | D                   | 0         | D        | 1                    | D                   | 3,675                  | H                     | -2,3         | D           | -7,86         | D            | 0,979       | 7,514    | 34         | 0,997      | 0,985                   | D                      | 0,946         | D            | 0,863        | D           | 0,707                    | 0                           | 5,76    | 0,791                 | 0,935                 | 0,998                    | 1                        | 18,523              | Deleterious            |
| 1  | chr3  | 128204951 | C   | T   | SNV      | Coding,Intergenic                          | missense_variant,upstream_gene_variant                                                              | rs2335052            | COSM445531  |                | 100  | 54,9         | 1312        | 587       | 720       | -100        | chr3  | 128204951 | 128204951 | C   | T   | exonic       | GATA2        | 0                  | nonsynonymous SNV   | GATA2:NM_001145662:exon3:c,G490A:p,A164T,GATA2:NM_032638:exon3:c,G490A:p,A164T,GATA2:NM_001145661:exon4:c,G490A:p,A164T                                                                                                                                                                                                                                | 3q21,3   | 0,1638          | rs2335052   | 0,429      | T         | 0,013                | B                   | 0,025                | B                   | 0,004     | N        | 0,729                | P                   | -0,175                 | N                     | -4,32        | D           | -0,43         | N            | 0,117       | 1,482    | 13,22      | 0,99       | 0,91                    | D                      | -1,065        | T            | 0            | T           | 0,543                    | 0                           | 1,83    | 0,598                 | 0,739                 | 0,999                    | 0,999                    | 9,555               | Tolerated              |
| 2  | chr3  | 128205860 | G   | C   | SNV      | Coding,Intergenic                          | synonymous_variant,upstream_gene_variant                                                            | rs1573858            |             |                | 100  | 99,7         | 1925        | 0         | 1920      | -100        | chr3  | 128205860 | 128205860 | G   | C   | exonic       | GATA2        | 0                  | synonymous SNV      | GATA2:NM_001145662:exon2:c,C15G:p,P5P,GATA2:NM_032638:exon2:c,C15G:p,P5P,GATA2:NM_001145661:exon3:c,C15G:p,P5P                                                                                                                                                                                                                                         | 3q21,3   | 0,6934          | rs1573858   | 0          | 0         | 0                    | 0                   | 0                    | 0                   | 0         | 0        | 0                    | 0                   | 0                      | 0                     | 0            | 0           | 0             | 0            | 0           | 0        | 0          | 0          | 0                       | 0                      | 0             | 0            | 0            | 0           | 0                        | 0                           | 0       | 0                     | 0                     | 0                        | 0                        | 0                   | Tolerated (synonymous) |
| 3  | chr4  | 55141055  | A   | G   | SNV      | Coding                                     | synonymous_variant                                                                                  | rs1873778            | COSM1430082 |                | 100  | 99,7         | 12713       | 35        | 12670     | -100        | chr4  | 55141055  | 55141055  | A   | G   | exonic       | PDGFRA       | 0                  | synonymous SNV      | PDGFRA:NM_006206:exon12:c,A1701G:p,P567P                                                                                                                                                                                                                                                                                                               | 4q12     | 0,9589          | rs1873778   | 0          | 0         | 0                    | 0                   | 0                    | 0                   | 0         | 0        | 0                    | 0                   | 0                      | 0                     | 0            | 0           | 0             | 0            | 0           | 0        | 0          | 0          | 0                       | 0                      | 0             | 0            | 0            | 0           | 0                        | 0                           | 0       | 0                     | 0                     | 0                        | 0                        | 0                   | Tolerated (synonymous) |
| 4  | chr4  | 106196829 | T   | G   | SNV      | Coding                                     | missense_variant                                                                                    | rs34402524           | COSM87176   |                | 100  | 48,2         | 38535       | 19801     | 18585     | -100        | chr4  | 106196829 | 106196829 | T   | G   | exonic       | TET2         | 0                  | nonsynonymous SNV   | TET2:NM_001127208:exon11:c,T5162G:p,L1721W                                                                                                                                                                                                                                                                                                             | 4q24     | 0,1233          | rs34402524  | 0,003      | D         | 0,794                | P                   | 0,754                | P                   | 0         | 0        | 1                    | P                   | 0                      | N                     | 4,35         | T           | 0,03          | N            | 0,232       | 2,834    | 21,5       | 0,932      | 0,78                    | D                      | -0,626        | T            | 0            | T           | 0,707                    | 0                           | 5,16    | 0,991                 | 1,011                 | 0,775                    | 0,124                    | 13,59               | Tolerated              |
| 5  | chr4  | 106196937 | T   | -0  | Deletion | Coding                                     | frameshift_variant,feature_truncation                                                               |                      |             |                | 100  | 38,5         | 31380       | 19297     | 12083     | -100        | chr4  | 106196937 | 106196937 | T   | -0  | exonic       | TET2         | 0                  | frameshift deletion | TET2:NM_001127208:exon11:c,5270delA:p,H1757fs                                                                                                                                                                                                                                                                                                          | 4q24     | 0               | 0           | 0          | 0         | 0                    | 0                   | 0                    | 0                   | 0         | 0        | 0                    | 0                   | 0                      | 0                     | 0            | 0           | 0             | 0            | 0           | 0        | 0          | 0          | 0                       | 0                      | 0             | 0            | 0            | 0           | 0                        | 0                           | 0       | 0                     | 0                     | 0                        | 0                        | 0                   | 0                      |
| 6  | chr7  | 101917521 | G   | A   | SNV      |                                            |                                                                                                     | rs803064             |             |                | 100  | 55,8         | 9023        | 3973      | 5033      | -100        | chr7  | 101917521 | 101917521 | G   | A   | exonic       | CUX1         | 0                  | nonsynonymous SNV   | CUX1:NM_001202544:exon15:c,G1342A:p,A448T,CUX1:NM_001202545:exon15:c,G1252A:p,A418T,CUX1:NM_001202546:exon15:c,G1273A:p,A425T,CUX1:NM_001913:exon16:c,G1390A:p,A464T,CUX1:NM_181500:exon16:c,G1384A:p,A462T                                                                                                                                            | 7q22,1   | 0,575           | rs803064    | 0,902      | T         | 0,229                | B                   | 0,052                | B                   | 0         | 0        | 1                    | P                   | 1,24                   | L                     | 1,51         | T           | -1,17         | N            | 0,273       | 0,313    | 5,82       | 0,736      | 0,84                    | D                      | -0,953        | T            | 0            | T           | 0,706                    | 0                           | -0,106  | -0,683                | -0,839                | 0,437                    | 0                        | 7,499               | Tolerated              |
| 7  | chr7  | 101921289 | A   | G   | SNV      |                                            |                                                                                                     | rs2230103,rs79334660 |             |                | 100  | 37           | 2903        | 1830      | 1073      | -100        | chr7  | 101921289 | 101921289 | A   | G   | exonic       | CUX1         | 0                  | nonsynonymous SNV   | CUX1:NM_001202544:exon17:c,A1585G:p,I529V,CUX1:NM_001202545:exon17:c,A1495G:p,I499V,CUX1:NM_001202546:exon17:c,A1516G:p,I506V,CUX1:NM_001913:exon18:c,A1633G:p,I545V,CUX1:NM_181500:exon18:c,A1627G:p,I543V                                                                                                                                            | 7q22,1   | 0,0382          | rs2230103   | 1          | T         | 0,139                | B                   | 0,07                 | B                   | 0         | 0        | 0,995                | N                   | 0,5                    | N                     | 1,63         | T           | 0,03          | N            | 0,326       | 0,53     | 7,655      | 0,604      | 0,967                   | D                      | -1,033        | T            | 0,004        | T           | 0,706                    | 0                           | 3,4     | 0,697                 | 1,186                 | 0,993                    | 1                        | 10,046              | Tolerated              |
| 8  | chr9  | 21971127  | A   | G   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,synonymous_variant                    |                      |             |                | 100  | 17,2         | 157         | 120       | 27        | -100        | chr9  | 21971127  | 21971127  | A   | G   | exonic       | CDKN2A       | 0                  | nonsynonymous SNV   | CDKN2A:NM_058195:exon2:c,T274C:p,S92P                                                                                                                                                                                                                                                                                                                  | 9p21,3   | 0               | 0           | 1          | T         | 0                    | B                   | 0                    | B                   | 0,471     | N        | 1                    | D                   | 0                      | 0                     | -0,88        | T           | 3,69          | N            | 0           | -0,153   | 1,37       | 0,714      | 0,169                   | N                      | -0,976        | T            | 0,192        | T           | 0,677                    | 0                           | 0,609   | -0,056                | 0,149                 | 0,536                    | 0,933                    | 4,699               | Tolerated              |
| 9  | chr9  | 21971127  | A   | T   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,synonymous_variant                    |                      |             |                | 44   | 5,7          | 157         | 120       | 9         | -100        | chr9  | 21971127  | 21971127  | A   | T   | exonic       | CDKN2A       | 0                  | nonsynonymous SNV   | CDKN2A:NM_058195:exon2:c,T274A:p,S92T                                                                                                                                                                                                                                                                                                                  | 9p21,3   | 0               | 0           | 0,236      | T         | 0,167                | B                   | 0,045                | B                   | 0,471     | N        | 1                    | D                   | 0                      | 0                     | -0,95        | T           | -0,45         | N            | 0           | -0,609   | 0,119      | 0,96       | 0,799                   | D                      | -0,99         | T            | 0,178        | T           | 0,677                    | 0                           | 0,609   | -0,056                | 0,149                 | 0,536                    | 0,933                    | 4,699               | Tolerated              |
| 10 | chr9  | 21971129  | T   | G   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,missense_variant                      |                      |             |                | 44   | 5,7          | 158         | 143       | 9         | -100        | chr9  | 21971129  | 21971129  | T   | G   | exonic       | CDKN2A       | 0                  | nonsynonymous SNV   | CDKN2A:NM_000077:exon2:c,A229C:p,T77P,CDKN2A:NM_001195132:exon2:c,A229C:p,T77P,CDKN2A:NM_058195:exon2:c,A272C:p,H91P                                                                                                                                                                                                                                   | 9p21,3   | 0               | 0           | 0,008      | D         | 0,998                | D                   | 0,915                | D                   | 0,029     | N        | 0,989                | N                   | 0,975                  | L                     | -0,99        | T           | -6,87         | D            | 0           | 4,105    | 23,7       | 0,982      | 0,938                   | D                      | -0,261        | T            | 0,379        | T           | 0,677                    | 0                           | 5,79    | 0,991                 | 1,061                 | 0,464                    | 0,913                    | 10,099              | Tolerated              |
| 11 | chr9  | 21971130  | G   | A   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,synonymous_variant                    |                      |             |                | 68   | 7,4          | 163         | 144       | 12        | -100        | chr9  | 21971130  | 21971130  | G   | A   | exonic       | CDKN2A       | 0                  | nonsynonymous SNV   | CDKN2A:NM_058195:exon2:c,C271T:p,H91Y                                                                                                                                                                                                                                                                                                                  | 9p21,3   | 0               | 0           | 0,197      | T         | 0,417                | B                   | 0,158                | B                   | 0,029     | N        | 1                    | N                   | 0                      | 0                     | -0,89        | T           | -3,72         | D            | 0           | 0,976    | 10,53      | 0,987      | 0,797                   | D                      | -0,845        | T            | 0,224        | T           | 0,677                    | 0                           | 3,74    | -0,711                | -0,197                | 0,328                    | 0,866                    | 5,67                | Tolerated              |
| 13 | chr17 | 7579472   | G   | C   | SNV      | Coding                                     | missense_variant                                                                                    | rs1042522            | COSM45985   | non-pathogenic | 100  | 59,1         | 3567        | 1452      | 2108      | -100        | chr17 | 7579472   | 7579472   | G   | C   | exonic       | TP53         | 0                  | nonsynonymous SNV   | TP53:NM_001126118:exon3:c,C98G:p,P33R,TP53:NM_000546:exon4:c,C215G:p,P72R,TP53:NM_001126112:exon4:c,C215G:p,P72R,TP53:NM_001126113:exon4:c,C215G:p,P72R,TP53:NM_001126114:exon4:c,C215G:p,P72R,TP53:NM_001276695:exon4:c,C98G:p,P33R,TP53:NM_001276696:exon4:c,C98G:p,P33R,TP53:NM_001276760:exon4:c,C98G:p,P33R,TP53:NM_001276761:exon4:c,C98G:p,P33R | 17p13,1  | 0,63            | rs1042522   | 0,642      | T         | 0,745                | P                   | 0,372                | B                   | 0,371     | U        | 1                    | P                   | 0                      | N                     | -2,05        | D           | -0,19         | N            | 0,267       | -0,415   | 0,355      | 0,57       | 0,361                   | N                      | -0,929        | T            | 0            | T           | 0,722                    | 0                           | 1,87    | 0,518                 | 1,045                 | 0,001                    | 0,001                    | 9,773               | Tolerated              |
| 12 | chr9  | 21971131  | G   | A   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,missense_variant                      |                      | COSM13766   |                | 44   | 5,7          | 157         | 144       | 9         | -100        | chr9  | 21971131  | 21971131  | G   | A   | exonic       | CDKN2A       | 0                  | nonsynonymous SNV   | CDKN2A:NM_000077:exon2:c,C227T:p,A76V,CDKN2A:NM_001195132:exon2:c,C227T:p,A76V                                                                                                                                                                                                                                                                         | 9p21,3   | 0               | 0           | 0,276      | T         | 0,004                | B                   | 0,001                | B                   | 0         | 0        | 1                    | N                   | 0                      | 0                     | -0,16        | T           | -1,74         | N            | 0,46        | 0,913    | 10,15      | 0,974      | 0,098                   | N                      | -1,03         | T            | 0,109        | T           | 0,677                    | 0                           | -0,271  | -0,011                | -0,239                | 0,309                    | 0,856                    | 6,3                 | Tolerated              |
| 14 | chr17 | 74733099  | G   | A   | SNV      | Intergenic,Intergenic,Coding,Intron,5P_UTR | downstream_gene_variant,upstream_gene_variant,synonymous_variant,intron_variant,5_prime_UTR_variant | rs237057             |             |                | 100  | 99,5         | 7558        | 34        | 7518      | -100        | chr17 | 74733099  | 74733099  | G   | A   | exonic       | SRSF2        | 0                  | synonymous SNV      | SRSF2:NM_001195427:exon1:c,C144T:p,D48D,SRSF2:NM_003016:exon1:c,C144T:p,D48D                                                                                                                                                                                                                                                                           | 17q25,1  | 0,797           | rs237057    | 0          | 0         | 0                    | 0                   | 0                    | 0                   | 0         | 0        | 0                    | 0                   | 0                      | 0                     | 0            | 0           | 0             | 0            | 0           | 0        | 0          | 0          | 0                       | 0                      | 0             | 0            | 0            | 0           | 0                        | 0                           | 0       | 0                     | 0                     | 0                        | 0                        | 0                   | Tolerated (synonymous) |
| 15 | chr20 | 31022959  | T   | C   | SNV      | Coding                                     | missense_variant                                                                                    | rs6058694            |             |                | 100  | 99,7         | 17877       | 39        | 17831     | -100        | chr20 | 31022959  | 31022959  | T   | C   | exonic       | ASXL1        | 0                  | nonsynonymous SNV   | ASXL1:NM_015338:exon12:c,T2444C:p,L815P                                                                                                                                                                                                                                                                                                                | 20q11,21 | 0,9999          | rs6058694   | 0          | 0         | 0                    | 0                   | 0                    | 0                   | 0         | 0        | 0                    | 0                   | 0                      | 0                     | 0            | 0           | 0             | 0            | 0           | 0        | 0          | 0          | 0                       | 0                      | 0             | 0            | 0            | 0           | 0                        | 0                           | 0       | 0                     | 0                     | 0                        | 0                        | 0                   | 0                      |
| 16 | chr21 | 36252877  | C   | T   | SNV      | Coding                                     | missense_variant                                                                                    |                      | COSM96546   |                | 100  | 46,5         | 8562        | 4573      | 3978      | -100        | chr21 | 36252877  | 36252877  | C   | T   | exonic       | RUNX1        | 0                  | nonsynonymous SNV   | RUNX1:NM_001001890:exon2:c,G404A:p,R135K,RUNX1:NM_001122607:exon2:c,G404A:p,R135K,RUNX1:NM_001754:exon5:c,G485A:p,R162K                                                                                                                                                                                                                                | 21q22,12 | 0               | 0           | 0,001      | D         | 0,999                | D                   | 0,997                | D                   | 0         | D        | 1                    | D                   | 3,195                  | M                     | -6,37        | D           | -2,5          | D            | 0,954       | 7,021    | 33         | 0,998      | 0,989                   | D                      | 0,994         | D            | 0,991        | D           | 0,722                    | 0                           | 5,31    | 0,871                 | 0,932                 | 1                        | 1                        | 19,335              | Deleterious            |
| 17 | chrX  | 39933339  | A   | G   | SNV      | Coding                                     | synonymous_variant                                                                                  | rs5917933            |             |                | 100  | 99,8         | 14616       | 18        | 14592     | -100        | chrX  | 39933339  | 39933339  | A   | G   | exonic       | BCOR         | 0                  | synonymous SNV      | BCOR:NM_001123383:exon4:c,T1260C:p,D420D,BCOR:NM_001123384:exon4:c,T1260C:p,D420D,BCOR:NM_001123385:exon4:c,T1260C:p,D420D,BCOR:NM_017745:exon4:c,T1260C:p,D420D                                                                                                                                                                                       | Xp11,4   | 0,8957          | rs5917933   | 0          | 0         | 0                    | 0                   | 0                    | 0                   | 0         | 0        | 0                    | 0                   | 0                      | 0                     | 0            | 0           | 0             | 0            | 0           | 0        | 0          | 0          | 0                       | 0                      | 0             | 0            | 0            | 0           | 0                        | 0                           | 0       | 0                     | 0                     | 0                        | 0                        | 0                   | Tolerated (synonymous) |
| 18 | chrX  | 123195650 | A   | T   | SNV      | Coding                                     | missense_variant                                                                                    |                      |             |                | 47   | 7,3          | 109         | 100       | 8         | -100        | chrX  | 123195650 | 123195650 | A   | T   | exonic       | STAG2        | 0                  | nonsynonymous SNV   | STAG2:NM_006603:exon16:c,A1564T:p,I522F,STAG2:NM_001042749:exon17:c,A1564T:p,I522F,STAG2:NM_001042750:exon17:c,A1564T:p,I522F,STAG2:NM_001042751:exon17:c,A1564T:p,I522F,STAG2:NM_001282418:exon17:c,A1564T:p,I522F                                                                                                                                    | Xq25     | 0               | 0           | 0,001      | D         | 0,955                | P                   | 0,774                | P                   | 0         | D        | 1                    | D                   | 2,66                   | M                     | 1,51         | T           | -3,81         | D            | 0,906       | 6,183    | 28,6       | 0,989      | 0,973                   | D                      | -0,343        | T            | 0,303        | T           | 0                        | 0                           | 4,77    | 1,062                 | 1,088                 | 0,999                    | 1                        | 13,668              | Tolerated              |

Known bugs and missing features
-------------
Some deletions are not handled correctly and sometimes there is no score for the deletion available.
Also there is no scoring system for non coding variants, this may change later. 

References
-------------
<sup>1</sup> http://annovar.openbioinformatics.org/en/latest/
https://www.gesundheit.gv.at/labor/laborwerte/blutbild/blasten
http://annovar.openbioinformatics.org/
http://cadd.gs.washington.edu/download
http://www.ensembl.org/info/genome/variation/predicted_data.html
http://varianttools.sourceforge.net/Annotation/DbNSFP
https://brb.nci.nih.gov/seqtools/colexpanno.html#dbnsfp
http://www.illumina.com/products/by-type/clinical-research-products/trusight-myeloid.html
http://www.enlis.com/blog/2015/03/17/the-best-variant-prediction-method-that-no-one-is-using/

Literatur
-----------
- Papaemmanuil, Elli; Gerstung, Moritz; Bullinger, Lars; Gaidzik, Verena I.; Paschka, Peter; Roberts, Nicola D. et al. (2016): Genomic Classification and Prognosis in Acute Myeloid Leukemia. In: The New England journal of medicine 374 (23), S. 2209–2221. DOI: 10.1056/NEJMoa1516192         
- Yang, Hui; Wang, Kai (2015): Genomic variant annotation and prioritization with ANNOVAR and wANNOVAR. In: Nature protocols 10 (10), S. 1556–1566. DOI: 10.1038/nprot.2015.105          
- Liu, Xiaoming; Jian, Xueqiu; Boerwinkle, Eric (2011): dbNSFP: a lightweight database of human nonsynonymous SNPs and their functional predictions. In: Human mutation 32 (8), S. 894–899. DOI: 10.1002/humu.21517          
- Dong, Chengliang; Wei, Peng; Jian, Xueqiu; Gibbs, Richard; Boerwinkle, Eric; Wang, Kai; Liu, Xiaoming (2015): Comparison and integration of deleteriousness prediction methods for nonsynonymous SNVs in whole exome sequencing studies. In: Human molecular genetics 24 (8), S. 2125–2137. DOI: 10.1093/hmg/ddu733   
- Wakita, S.; Yamaguchi, H.; Ueki, T.; Usuki, K.; Kurosawa, S.; Kobayashi, Y. et al. (2016): Complex molecular genetic abnormalities involving three or more genetic mutations are important prognostic factors for acute myeloid leukemia. In: Leukemia 30 (3), S. 545–554. DOI: 10.1038/leu.2015.288.
   



> Written with [StackEdit](https://stackedit.io/).
