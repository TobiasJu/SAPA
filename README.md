SAPA - SNP annotation programm for AML
===================

This program is designed to add additional information to an Illumina truseq amplicon variants file.

The following Parameters are added to the file:

- COMING SOON




#### <i class="icon-file"></i> Input File
As input file is a Illumina truseq amplicon variants file needed, which contains the patients SNPs in a character separated format.


Please Note that during the first run of the program, the required database will be downloaded. This may take some time, depending on your internet connection.

Several commonly used databases are integrated: ‘cytoBand’ for the chromosome coordinate of each cytogenetic band,
‘1000g2014oct’ for alternative allele frequency in the 1000 Genomes Project (version October 2014), ‘exac03’
for the variants reported in the Exome Aggregation Consortium (version 0.3)50, ‘ljb26_all’ for various functional
deleteriousness prediction scores from the dbNSFP database (version 2.6)51, ‘clinvar_20140929’ for the variants
reported in the ClinVar database (version 20140929)52 and ‘snp138’ for the dbSNP database (version 138)53.
