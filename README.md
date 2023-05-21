# TFG - Methylation analysis of melon varieties using nanopore sequencing
Code used in the method and results obtained for the analysis of methylation, annotation of unannotated genes of _Cucumis melo_ from _Arabidopsis thaliana_ and _Vitis vinifera_ and relationship of total annotated genes with processes of synthesis of volatile organic compounds (VOCs) developed in the Final Degree Project at the University of Murcia during the academic year 2022/2023.
## Scripts
-  **Analysis.R**: From the results obtained from the nanoporous sequencing in TSV format, data manipulation is performed to obtain the methylation study of the different varieties separately and as a whole, allowing the analysis of coincident and unique methylations of each variety among other characteristics. 
-  **tBlastn_Annotation.R**: From the unannotated genes of each variety of _Cucumis melo_, its annotations are retrieved by querying Biomart and running tBlastn. In addition, from the annotations retrieved, it selects which of them are chosen for each gene. It also contains the necessary functionality to measure the performance of the method through the calculation of different metrics.
-  **BP_Research.R**: After obtaining the annotations of the unannotated methylated genes, the previously annotated genes are taken together with the genes whose annotations have been found using the bioinformatics method and a search is performed for those GO annotations that share with genes involved in the synthesis of VOCs.
## Dependencies
Scripts use one or more of the following libraries:
- **Tidyverse**: Set of R packages designed for data science. It allows importing, transforming, visualizing, modeling and communicating data sets.
- **biomaRt**: Package that provides an interface to a collection of databases that implement the BioMart software package. In this case, Ensembl Plants.
- **openxlsx**: Package that allows you to read, write and edit xlsx files.
- **seqinr**: Package for data analysis and visualization.
- **dplyr**:  Package that provides a "grammar" for manipulation and operations with data frames.
## Starting files
The starting files are those obtained from nanopore sequencing. There is one for each variety:
- **methylation_info_AM3-6.tsv**
- **methylation_info_CH2.tsv**
- **methylation_info_PSA2.tsv**
- **methylation_info_PSU2.tsv**
- **methylation_info_ROCH2.tsv**
- **methylation_info_TN2.tsv**
- **methylation_info_TB2.tsv**
## FASTA files of genomes
The FASTA files containing the genome of the different organisms used for the execution of tBlastn are also required..
- **_Arabidopsis thaliana_**: https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_blastsets/TAIR10_cdna_20101214_updated
- _**Solanum lycopersicum**_: https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_CDS.fasta
- _**Vitis vinifera**_: https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/fasta/vitis_vinifera/cds/Vitis_vinifera.PN40024.v4.cds.all.fa.gz
- _**Amborella trichopoda**_: https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/fasta/amborella_trichopoda/cdna/Amborella_trichopoda.AMTR1.0.cdna.all.fa.gz
