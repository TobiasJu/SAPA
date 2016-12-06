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

`-id, --input_directory` 
: hg19 (GRCh37) database directory (default /hg19)
	> e.g. `-id humandb/`
	set the input directory to be named humandb 

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
	- Score	ranges from 0 to 1. The amino acid substitution is predicted damaging is the score is <= 0.05, and tolerated if the score is > 0.05.
	- http://sift.jcvi.org/www/SIFT_help.html
	
- **Polyphen2 **
	- Scores ranges from 0 to 1
	- D: Probably damaging (>=0.957), P: possibly damaging (0.453<=pp2_hdiv<=0.956); B: benign (pp2_hdiv<=0.452)
	- http://genetics.bwh.harvard.edu/pph2/dokuwiki/overview

- **LRT ** 
	- Scores ranges from 0 to 1
	- D: Deleterious; N: Neutral; U: Unknown
	- https://www.researchgate.net/profile/Justin_Fay/publication/26671245_Identification_of_deleterious_mutations_within_three_human_genomes/links/00b49529eabbc6aa2b000000.pdf

- **MutationTaster **
	- Scores ranges from 0 to 1
	- 	A = "disease_causing_automatic"; D =“disease_causing”; N = “polymorphism”; P = “polymorphism_automatic”
	- http://www.mutationtaster.org/


- **MutationAssessor**
	- Scores ranges from -5.545 to 5.975
	- H: high; M: medium; L: low; N: neutral. H/M means functional and L/N means non-functional
	- http://mutationassessor.org/r3/


- **FATHMM**
	- Scores ranges from -16.13 to 10.64 
	- D: Deleterious; T: Tolerated
	- http://fathmm.biocompute.org.uk/


- **PROVEAN **
	- Scores ranges from -13 to 4
	- Deleterious PROVEAN < -2,5 ; Tolerated PROVEAN > -2,5 
	- http://provean.jcvi.org/about.php


- **VEST3**
	- Scores ranges from 0 to 1. 
	- The larger the score the more likely the mutation may cause functional change.
	- http://karchinlab.org/apps/appVest.html


- **CADD**
	- Scores ranges from ?? to ??
	- The larger the score the more likely the SNP has damaging effect.
	- http://cadd.gs.washington.edu/

- **CADD_phred**
	- This is phred-like rank score based on whole genome CADD raw scores.
	- The larger the score the more likely the SNP has damaging effect.

- **DANN**
	- Scores ranges from 0.01 to 0.99
	- The higher the score the more likely the mutation may cause functional change.
	- DANN uses the same feature set and training data as CADD to train a deep neural network (DNN)
	- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4341060/

- **FATHMM_MKL_coding_score**
	- Scores ranges from 0 - 1
	- The higher the score, the greater the confidence of a functional mutation.
	- http://fathmm.biocompute.org.uk/fathmmMKL.htm

- **MetaSVM**
	- Scores ranges from -1 to +1
	- D: Deleterious; T: Tolerated
	- Uses a radial SVM model to train prediction model, using all available scoring algorithm normalized scores
	- http://annovar.openbioinformatics.org/en/latest/user-guide/filter/#-metasvm-annotation

- **MetaLR **
	- Scores ranges from 0 to 1
	- D: Deleterious; T: Tolerated
	- very similar to MetaSVM, with similar performance
	- http://annovar.openbioinformatics.org/en/latest/user-guide/filter/#-metalr-annotation

- **integrated fitCons**
	- Scores ranges from 0.0 to 0.80
	- http://compgen.cshl.edu/fitCons/

- **GERP RS**
	- Scores ranges from -12.3 to 6.17
	- The larger the score, the more conserved the site.
	- http://mendel.stanford.edu/SidowLab/downloads/gerp/

- **phyloP7way_vertebrate**
	- The larger the score, the more conserved.
		- higher scores are more deleterious
	- Phylogenetic p-values for 7 vertebrate species
	- https://brb.nci.nih.gov/seqtools/colexpanno.html#dbnsfp

- **phyloP20way_mammalian**
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

#### <i class="icon-file"></i> Output File (default)
The generated outputfile `output.txt`, which contains the patients annotated SNPs in a character separated format, will be in the current working directory. Using the included example dataset in the /data directory, the following output will be generated.

| ID | Chr   | Pos       | Ref | Alt | Type     | Context                                    | Consequence                                                                                         | dbSNP                | COSMIC      | ClinVar        | Qual | Alt Freq | Total Depth | ,Ref Depth | Alt Depth | Strand Bias | function prediction scores | conservation scores | ensemble scores | final prediction |
|----|-------|-----------|-----|-----|----------|--------------------------------------------|-----------------------------------------------------------------------------------------------------|----------------------|-------------|----------------|------|----------|-------------|------------|-----------|-------------|----------------------------|---------------------|-----------------|------------------|
| 0  | chr2  | 25463566  | C   | A   | SNV      | Coding                                     | missense_variant                                                                                    |                      |             |                | 100  | 334      | 6568        | 4358       | 2194      | -100.0      | 863                        | 5.76                | 4810            | Deleterious      |
| 1  | chr3  | 128204951 | C   | T   | SNV      | Coding,Intergenic                          | missense_variant,upstream_gene_variant                                                              | rs2335052            | COSM445531  |                | 100  | 549      | 1312        | 587        | 720       | -100.0      | 0                          | 1.83                | 1801            | Tolerated        |
| 2  | chr3  | 128205860 | G   | C   | SNV      | Coding,Intergenic                          | synonymous_variant,upstream_gene_variant                                                            | rs1573858            |             |                | 100  | 997      | 1925        | 0          | 1920      | -100.0      | .                          | .                   | .               | .                |
| 3  | chr4  | 55141055  | A   | G   | SNV      | Coding                                     | synonymous_variant                                                                                  | rs1873778            | COSM1430082 |                | 100  | 997      | 12713       | 35         | 12670     | -100.0      | .                          | .                   | .               | .                |
| 4  | chr4  | 106196829 | T   | G   | SNV      | Coding                                     | missense_variant                                                                                    | rs34402524           | COSM87176   |                | 100  | 482      | 38535       | 19801      | 18585     | -100.0      | 0                          | 5.16                | 1344            | Tolerated        |
| 5  | chr4  | 106196937 | T   | -0  | Deletion | Coding                                     | frameshift_variant,feature_truncation                                                               |                      |             |                | 100  | 385      | 31380       | 19297      | 12083     | -100.0      | .                          | .                   | .               | .                |
| 6  | chr7  | 101917521 | G   | A   | SNV      |                                            |                                                                                                     | rs803064             |             |                | 100  | 558      | 9023        | 3973       | 5033      | -100.0      | 0                          | -106                | 489             | Tolerated        |
| 7  | chr7  | 101921289 | A   | G   | SNV      |                                            |                                                                                                     | rs2230103,rs79334660 |             |                | 100  | 0.37     | 2903        | 1830       | 1073      | -100.0      | 4                          | 3.4                 | 90              | Tolerated        |
| 8  | chr9  | 21971127  | A   | G   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,synonymous_variant                    |                      |             |                | 100  | 172      | 157         | 120        | 27        | -100.0      | 192                        | 609                 | 620             | Tolerated        |
| 9  | chr9  | 21971127  | A   | T   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,synonymous_variant                    |                      |             |                | 44   | 57       | 157         | 120        | 9         | -100.0      | 178                        | 609                 | 2835            | Tolerated        |
| 10 | chr9  | 21971129  | T   | G   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,missense_variant                      |                      |             |                | 44   | 57       | 158         | 143        | 9         | -100.0      | 379                        | 5.79                | 3201            | Tolerated        |
| 11 | chr9  | 21971130  | G   | A   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,synonymous_variant                    |                      |             |                | 68   | 74       | 163         | 144        | 12        | -100.0      | 224                        | 3.74                | 2090            | Tolerated        |
| 13 | chr17 | 7579472   | G   | C   | SNV      | Coding                                     | missense_variant                                                                                    | rs1042522            | COSM45985   | non-pathogenic | 100  | 591      | 3567        | 1452       | 2108      | -100.0      | 0                          | 1.87                | 823             | Tolerated        |
| 12 | chr9  | 21971131  | G   | A   | SNV      | Intron,Intergenic,Coding                   | intron_variant,NMD_transcript_variant,downstream_gene_variant,missense_variant                      |                      | COSM13766   |                | 44   | 57       | 157         | 144        | 9         | -100.0      | 109                        | -271                | 1522            | Tolerated        |
| 14 | chr17 | 74733099  | G   | A   | SNV      | Intergenic,Intergenic,Coding,Intron,5P_UTR | downstream_gene_variant,upstream_gene_variant,synonymous_variant,intron_variant,5_prime_UTR_variant | rs237057             |             |                | 100  | 995      | 7558        | 34         | 7518      | -100.0      | .                          | .                   | .               | .                |
| 15 | chr20 | 31022959  | T   | C   | SNV      | Coding                                     | missense_variant                                                                                    | rs6058694            |             |                | 100  | 997      | 17877       | 39         | 17831     | -100.0      | 0                          | -4.18               | 652             | Tolerated        |
| 16 | chr21 | 36252877  | C   | T   | SNV      | Coding                                     | missense_variant                                                                                    |                      | COSM96546   |                | 100  | 465      | 8562        | 4573       | 3978      | -100.0      | 991                        | 5.31                | 5748            | Deleterious      |
| 17 | chrX  | 39933339  | A   | G   | SNV      | Coding                                     | synonymous_variant                                                                                  | rs5917933            |             |                | 100  | 998      | 14616       | 18         | 14592     | -100.0      | .                          | .                   | .               | .                |
| 18 | chrX  | 123195650 | A   | T   | SNV      | Coding                                     | missense_variant                                                                                    |                      |             |                | 47   | 73      
 | 109         | 100        | 8         | -100.0      | 303                        | 4.77                | 4932            | Tolerated        |


Known bugs and missing features
-------------
Some deletions are not handled correctly and sometimes there is no score for the deletion available.
Also there is no scoring system for non coding variants, this may change later. 

References
-------------
<sup>1</sup> http://annovar.openbioinformatics.org/en/latest/

> Written with [StackEdit](https://stackedit.io/).
