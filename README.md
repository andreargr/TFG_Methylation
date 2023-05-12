# TFG - Análisis de metilación de variedades de melón usando secuenciación por nanoporos
Código utilizado en el método y obtención de resultados para el análisis de metilación, anotación de genes no anotados de _Cucumis melo_ a partir de _Arabidopsis thaliana_ y _Vitis vinifera_ y relación de genes anotados totales con procesos de síntesis de compuestos orgánicos volátiles (VOCs) desarrollado en el Trabajo de Fin de Grado en la Universidad de Murcia durante el curso 2022/2023.
## Scripts
-  **Analysis.R**: A partir de los resultados obtenidos de la secuenciación de nanoporos en formato TSV, se realiza la manipulación de datos para obtener el estudio de metilación de las distintas variedades por separado y en conjunto, permitiendo analizar metilaciones coincidentes y únicas de cada variedad entre otras características. 
-  **tBlastn_Annotation.R**: A partir de los genes no anotados de cada variedad de _Cucumis melo_, se obtienen sus anotaciones mediante consultas a Biomart y realización de tBlastn. Además, de las anotaciones que se recuperan, selecciona cuáles de ellas se escogen para cada gen. También, contiene la funcionalidad necesaria para medir el rendimiento del método a través del cálculo de distintas métricas.
-  **BP_Research.R**: Una vez se tienen los genes con metilaciones únicas anotados, busca si entre las distintas anotaciones se encuentran aquellas relacionadas con procesos de síntesis de VOCs.
## Dependencias
Los scripts utilizan alguna o varias de las siguientes librerías:
- **Tidyverse**: Conjunto de paquetes en R diseñados para ciencia de datos. Permite importar, transformar visualizar, modelar y comunicar conjuntos de datos.
- **biomaRt**: Paquete que proporciona una interfaz a una colección de bases de datos que implementan el paquete de software BioMart (). En este caso, Ensembl Plants.
- **openxlsx**: Paquete que permite leer, escribir y editar archivos xlsx.
- **seqinr**: Paquete destinado al análisis y visualización de datos.
- **dplyr**:  Paquete que proporciona una "gramática" para la manipulación y operaciones con data frames.
## Ficheros de partida
Los ficheros de partida son aquellos obtenidos a partir de la secuenciación de nanoporos. Hay uno para cada variedad:
- **methylation_info_AM3-6.tsv**
- **methylation_info_CH2.tsv**
- **methylation_info_PSA2.tsv**
- **methylation_info_PSU2.tsv**
- **methylation_info_ROCH2.tsv**
- **methylation_info_TN2.tsv**
- **methylation_info_TB2.tsv**
## Ficheros FASTA de los genomas
También son necesarios los ficheros FASTA que contienen el genoma de los distintos organismos utilizados para  la ejecución del tBlastn.
- **_Arabidopsis thaliana_**: https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_blastsets/TAIR10_cdna_20101214_updated
- _**Solanum lycopersicum**_: https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_CDS.fasta
- _**Vitis vinifera**_: https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/fasta/vitis_vinifera/cds/Vitis_vinifera.PN40024.v4.cds.all.fa.gz
- _**Amborella trichopoda**_: https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/fasta/amborella_trichopoda/cdna/Amborella_trichopoda.AMTR1.0.cdna.all.fa.gz
