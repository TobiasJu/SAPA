SAPA - SNP annotation programm for AML
===================

**75%** of a **Acute myeloid leukemia** (AML) patients SNPs are unique to them and their impact is still unresolved. 

SAPA is designed to add additional information to an **Illumina truseq amplicon** variants csv file. Such as Scores for nonsynonymous scoring matrices, divided into **3** different **combined Scores**, function prediction scores, conservation scores and one ensemble score. So that the viewer gets more Information whether a mutation is delirious or harmless. Due to performance and data reasons, SAPA uses annovar<sup>1</sup> and their precomputed SNP database, to quickly gather the deleterious prediction methods scores. 

[^annovar]: Annovar is available here: http://annovar.openbioinformatics.org/. 
  

First run 
-------------
<i class="icon-refresh"></i> To get startet simply type in:

`./main.py -i data/truseq_example_data.csv`

this will run the program with the example dataset. Your output will be saved in the `output.txt ` file.  
  

> **Please Note:** During the first run of the program, the required **database will be downloaded**. This may take some time, depending on your internet connection. You may want to use the `-fast` argument to minimize the download size.


Parameter
-------------
<a href="" rel="this is just cause of poor markdown implementation"><img src="images/paramters.png" alt="" /></a>
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


The following scores are added to the advanced output file:

- SIFT
- Polyphen2_HDIV
- LRT	MutationTaster	
- MutationAssessor
- FATHMM
- Radial SVM
- LR
- VEST3
- CADD_raw
- CADD_phred
- GERP_RS
- phyloP46way placental
- phyloP100way_vertebrate	
- SiPhy_29way_logOdds

For a complete description of the single scores visit:
http://annovar.openbioinformatics.org/en/latest/user-guide/filter/#ljb42-dbnsfp-non-synonymous-variants-annotation

#### <i class="icon-hdd"></i> Databases
All required databases will be downloaded in the hg19/ <i class="icon-folder-open"></i> directory. 
Several commonly used databases are integrated: ‘cytoBand’ for the chromosome coordinate of each cytogenetic band, ‘exac03’ for the variants reported in the Exome Aggregation Consortium (version 0.3), ‘ljb26_all’ for various functional deleteriousness prediction scores from the dbNSFP database (version 2.6), ‘clinvar_20140929’ for the variants reported in the ClinVar database (version 20140929).


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

References
-------------
<sup>1</sup> http://annovar.openbioinformatics.org/en/latest/

> Written with [StackEdit](https://stackedit.io/).