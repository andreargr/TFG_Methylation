library(dplyr)
library(seqinr)
library(biomaRt)

setwd("C:/Users/Andrea/Desktop/CUARTO/TFG/metilacion/basedatos")

ensembl_plants <- useEnsemblGenomes(biomart = "plants_mart")

#Anotaciones vitis
searchDatasets(ensembl_plants, pattern = "Vitis vinifera")
ensembl_vitis <- useEnsemblGenomes(biomart = "plants_mart", 
                                   dataset = "vvinifera_eg_gene")
#Anotaciones arabidopsis
searchDatasets(ensembl_plants, pattern = "thaliana")
ensembl_arabidopsis <- useEnsemblGenomes(biomart = "plants_mart", 
                                         dataset = "athaliana_eg_gene")
#tblastn
busca_blast <- function (query, bd, salida){
  query_string <- str_c("tblastn -query ", query, " -db ", bd, " -outfmt \"7 qacc sacc pident evalue qcovs qstart qend sstart send\"  -out ", salida)
  blast_out <- system(command = query_string)
  
}

#makeblastdb -in TAIR10_cdna_20101214_updated.fasta -dbtype nucl -out tair
#makeblastdb -in ITAG4.1_CDS.fasta -dbtype nucl -out tomato
#makeblastdb -in Vitis_vinifera.PN40024.v4.cds.all.fa -dbtype nucl -out uva
#makeblastdb -in Amborella_trichopoda.AMTR1.0.cdna.all.fa -dbtype nucl -out Amborella_tri

#Selección de los 75 genes de prueba
gnoanotados_cucumis <- gnoanotados_cucumis %>% distinct(Nearest.Unigene.x, .keep_all = TRUE) %>% dplyr::select("Nearest.Unigene.x","peptide") 
gnoanotados_cucumis <- gnoanotados_cucumis[1:50,]
write.fasta(sequences = as.list(gnoanotados_cucumis$peptide), names = gnoanotados_cucumis$Nearest.Unigene.x,file.out = "50genes.fasta")

tomate <- busca_blast("75genes.fasta","tomato", "tomate.txt")
tair <- busca_blast("75genes.fasta","tair", "tair.txt")
uva <- busca_blast("75genes.fasta","uva", "uva.txt")
amborella <- busca_blast("75genes.fasta","Amborella_tri", "amborella.txt")

nombres_columna <- c("qacc", "sacc","pidentity", "evalue","qcovs", "qstart", "qend", "sstart", "send")

#Filtro por evalue y %identidad
MAX_EV <- 0.000001
MIN_PID <- 70
MIN_PID2 <- 50
MIN_COV <- 80
MIN_COV2 <- 50
#CUCUMIS-SOLANUM LYCOPERSICUM
df_CMTomato <- read.table("tomate.txt",sep='\t', header = FALSE, col.names = nombres_columna)
filtrado_CMTomato <- df_CMTomato %>% filter(evalue <= MAX_EV) %>% filter(pidentity>=MIN_PID) %>% arrange(evalue)
filtrado_70id80covTomato <- df_CMTomato %>% filter(pidentity >= MIN_PID) %>% filter(qcovs>=MIN_COV) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
filtrado_70id50covTomato <- df_CMTomato %>% filter(pidentity >= MIN_PID) %>% filter(qcovs>=MIN_COV2) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
filtrado_50id80covTomato <- df_CMTomato %>% filter(pidentity >= MIN_PID2) %>% filter(qcovs>=MIN_COV) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
filtrado_50id50covTomato <- df_CMTomato %>% filter(pidentity >= MIN_PID2) %>% filter(qcovs>=MIN_COV2) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)

#CUCUMIS-ARABIDOPSIS THALIANA
df_CMArabidopsis <- read.table("tair.txt",sep='\t', header = FALSE, col.names = nombres_columna)
filtrado_CMArabidopsis <- df_CMArabidopsis %>% filter(evalue <= MAX_EV)%>% filter(pidentity>=MIN_PID)%>% arrange(evalue)

filtrado_70id80covArabidopsis <- df_CMArabidopsis %>% filter(pidentity >= MIN_PID) %>% filter(qcovs>=MIN_COV) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
filtrado_70id50covArabidopsis <- df_CMArabidopsis %>% filter(pidentity >= MIN_PID) %>% filter(qcovs>=MIN_COV2) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
filtrado_50id80covArabidopsis <- df_CMArabidopsis %>% filter(pidentity >= MIN_PID2) %>% filter(qcovs>=MIN_COV) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
filtrado_50id50covArabidopsis <- df_CMArabidopsis %>% filter(pidentity >= MIN_PID2) %>% filter(qcovs>=MIN_COV2) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
#CUCUMIS-VITIS VINIFERA
df_CMVitis <- read.table("uva.txt",sep='\t', header = FALSE, col.names = nombres_columna)
filtrado_CMVitis <- df_CMVitis %>% filter(evalue <= MAX_EV)%>% filter(pidentity>=MIN_PID)%>% arrange(evalue)

filtrado_70id80covVitis <- df_CMVitis %>% filter(pidentity >= MIN_PID) %>% filter(qcovs>=MIN_COV) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
filtrado_70id50covVitis <- df_CMVitis %>% filter(pidentity >= MIN_PID) %>% filter(qcovs>=MIN_COV2) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
filtrado_50id80covVitis <- df_CMVitis %>% filter(pidentity >= MIN_PID2) %>% filter(qcovs>=MIN_COV) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
filtrado_50id50covVitis <- df_CMVitis %>% filter(pidentity >= MIN_PID2) %>% filter(qcovs>=MIN_COV2) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
#CUCUMIS-AMBORELLA TRICHOPODA
df_CMAmborella <- read.table("amborella.txt",sep='\t', header = FALSE, col.names = nombres_columna)
filtrado_CMAmborella <- df_CMAmborella %>% filter(evalue <= MAX_EV)%>% filter(pidentity>=MIN_PID)%>% arrange(evalue)

filtrado_70id80covAmborella <- df_CMAmborella %>% filter(pidentity >= MIN_PID) %>% filter(qcovs>=MIN_COV) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
filtrado_70id50covAmborella <- df_CMAmborella %>% filter(pidentity >= MIN_PID) %>% filter(qcovs>=MIN_COV2) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
filtrado_50id80covAmborella <- df_CMAmborella %>% filter(pidentity >= MIN_PID2) %>% filter(qcovs>=MIN_COV) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)
filtrado_50id50covAmborella <- df_CMAmborella %>% filter(pidentity >= MIN_PID2) %>% filter(qcovs>=MIN_COV2) %>% filter(evalue <= MAX_EV) %>% arrange(evalue)

#graficas para determinar rango de coverage
MAX_EV <- 0.000001
#para todos los gnoanotados
tair <- busca_blast("gnoanotados.fasta","tair", "tair_all.txt") 
uva <- busca_blast("gnoanotados.fasta","uva", "uva_all.txt")
tomato <- busca_blast("gnoanotados.fasta","tomato", "tomato_all.txt")
amborella <- busca_blast("gnoanotados.fasta","Amborella_tri", "amborella_all.txt")

nombres_columna <- c("qacc", "sacc","pidentity", "evalue","qcovs", "qstart", "qend", "sstart", "send")

df_CMTomato_All <- read.table("tomato_all.txt",sep='\t', header = FALSE, col.names = nombres_columna)
df_CMArabidopsis_All <- read.table("tair_all.txt",sep='\t', header = FALSE, col.names = nombres_columna)
df_CMVitis_All <- read.table("uva_all.txt",sep='\t', header = FALSE, col.names = nombres_columna)
df_CMAmborella_All <- read.table("amborella_all.txt",sep='\t', header = FALSE, col.names = nombres_columna)

df_CMTomato_All2 <- df_CMTomato_All %>% arrange(evalue) %>% filter(evalue <= MAX_EV) %>% distinct(qacc, .keep_all = TRUE)
df_CMAmborella_All2 <- df_CMAmborella_All %>% arrange(evalue) %>% filter(evalue <= MAX_EV) %>% distinct(qacc, .keep_all = TRUE)
df_CMArabidopsis_All2 <- df_CMArabidopsis_All %>% arrange(evalue) %>% filter(evalue <= MAX_EV) %>% distinct(qacc, .keep_all = TRUE)
df_CMVitis_All2 <- df_CMVitis_All %>% arrange(evalue) %>% filter(evalue <= MAX_EV) %>% distinct(qacc, .keep_all = TRUE)

df_CMTomato_All2$especie <- "Tomato"
df_CMAmborella_All$especie <- "Amborella"
df_CMVitis_All2$especie <- "Vitis"
df_CMArabidopsis_All2$especie <- "Arabidopsis"

join_resultados <-  merge(df_CMTomato_All,df_CMAmborella_All, all=TRUE)
join_resultados <-  merge(join_resultados,df_CMArabidopsis_All, all=TRUE)
join_resultados_All <-  merge(join_resultados,df_CMVitis_All, all=TRUE)

gf_coverage_All <- join_resultados_All %>% ggplot(aes(x=especie,y=qcovs)) + geom_violin(colour="black", fill="yellow")
gf_id_All <-join_resultados_All %>% ggplot(aes(x=especie,y=pidentity)) + geom_violin(colour="black", fill="blue")


############################# ANOTACIÓN ########################
#Nos quedamos con Vitis vinifera y Arabidopsis thaliana, repetimos el proceso pero esta vez con todos los genes no anotados
gnoanotados_cucumis <- gnoanotados_cucumis %>% dplyr::select("Nearest.Unigene.x","peptide") 
write.fasta(sequences = as.list(gnoanotados_cucumis$peptide), names = gnoanotados_cucumis$Nearest.Unigene.x,file.out = "gnoanotados.fasta")

tair <- busca_blast("gnoanotados.fasta","tair", "tair_all.txt")
uva <- busca_blast("gnoanotados.fasta","uva", "uva_all.txt")

nombres_columna <- c("qacc", "sacc","pidentity", "evalue","qcovs", "qstart", "qend", "sstart", "send")

#####################################
#Para genes anotados - evaluación
ganotados_cucumis <- ganotados_cucumis %>% distinct(Nearest.Unigene.x, .keep_all = TRUE) #eliminamos repeticiones
write.fasta(sequences = as.list(ganotados_cucumis$peptide), names = ganotados_cucumis$Nearest.Unigene.x,file.out = "ganotados.fasta")

tair <- busca_blast("ganotados.fasta","tair", "tair_ganotados.txt")
uva <- busca_blast("ganotados.fasta","uva", "uva_ganotados.txt")

#####################################
#Filtro por evalue y %identidad
MAX_EV <- 0.000001
MIN_PID <- 20
#CUCUMIS-ARABIDOPSIS THALIANA genes no anotados
df_CMArabidopsis_all <- read.table("tair_all.txt",sep='\t', header = FALSE, col.names = nombres_columna)
filtrado_CMArabidopsis_all <- df_CMArabidopsis_all %>% filter(evalue <= MAX_EV)%>% filter(pidentity>=MIN_PID)%>% arrange(evalue)
filtrado_CMArabidopsis_all <- filtrado_CMArabidopsis_all %>% group_by(qacc) %>% distinct(qacc, .keep_all = TRUE) #nos quedamos con el resultado del gen con menor evalue

resultados_arabidopsis <- data.frame("ensembl_gene_id"=c(""), 'ensembl_transcript_id'=c(""), "go_id"=c(""),"description"=c(""),"peptide"=c(""))

for (i in 1:nrow(filtrado_CMArabidopsis_all)) {
  genCucumis <- filtrado_CMArabidopsis_all[i, "qacc"]
  genArabidopsis <- filtrado_CMArabidopsis_all[i, "sacc"]
  
  matches_arabidopsis <- list(genCucumis = list(genArabidopsis))
  for (genCucumis in matches_arabidopsis){
    genArabidopsis <- matches_arabidopsis[[1]][[1]]
    resultadosgenes_arabidopsis <- biomaRt::getBM(attributes = c("ensembl_gene_id",'ensembl_transcript_id',"go_id","description","peptide"),mart= ensembl_arabidopsis, filters = 'ensembl_transcript_id',values = unlist(genArabidopsis))
    resultados_arabidopsis <- full_join(resultados_arabidopsis,resultadosgenes_arabidopsis)
  }
}

resultados_arabidopsis <- resultados_arabidopsis[-c(1),] #eliminamos la fila vac?a
colnames(resultados_arabidopsis)[2]  <- "sacc" #cambiamos el nombre de la columna para realizar el leftjoin

resultados_arabidopsis <- resultados_arabidopsis %>% filter(go_id != c("NA","")) %>% left_join(filtrado_CMArabidopsis_all, by = c("sacc")) %>% dplyr::select(,-(9:12)) #eliminamos genes sin anotacion
resultados_arabidopsis <- resultados_arabidopsis %>% subset(go_id != "")
#####################################
#Para evaluación
ID50 <- 30
MAX_EV <- 0.000001
COV50 <- 50

df_CMArabidopsis_Evaluacion <- read.table("tair_ganotados.txt",sep='\t', header = FALSE, col.names = nombres_columna)
filtrado_Arabidopsis_Evaluacion <- df_CMArabidopsis_Evaluacion %>% arrange(evalue) %>% filter(evalue <= MAX_EV) %>% distinct(qacc, .keep_all = TRUE) %>% filter(pidentity>=ID50) %>% filter(qcovs>=COV50)

resultados_arabidopsis_Evaluacion <- data.frame("ensembl_gene_id"=c(""), 'ensembl_transcript_id'=c(""), "go_id"=c(""),"description"=c(""),"peptide"=c(""))

for (i in 1:nrow(filtrado_Arabidopsis_Evaluacion)) {
  genCucumis <- filtrado_Arabidopsis_Evaluacion[i, "qacc"]
  genArabidopsis <- filtrado_Arabidopsis_Evaluacion[i, "sacc"]
  
  matches_arabidopsis <- list(genCucumis = list(genArabidopsis))
  for (genCucumis in matches_arabidopsis){
    genArabidopsis <- matches_arabidopsis[[1]][[1]]
    resultadosgenes_arabidopsis <- biomaRt::getBM(attributes = c("ensembl_gene_id",'ensembl_transcript_id',"go_id","description","peptide"),mart= ensembl_arabidopsis, filters = 'ensembl_transcript_id',values = unlist(genArabidopsis))
    resultados_arabidopsis_Evaluacion <- full_join(resultados_arabidopsis_Evaluacion,resultadosgenes_arabidopsis)
  }
}

resultados_arabidopsis_Evaluacion <- resultados_arabidopsis_Evaluacion[-c(1),] #eliminamos la fila vac?a
colnames(resultados_arabidopsis_Evaluacion)[2]  <- "sacc" #cambiamos el nombre de la columna para realizar el leftjoin

resultados_arabidopsis_Evaluacion <- resultados_arabidopsis_Evaluacion %>% filter(go_id != c("NA","")) %>% left_join(filtrado_Arabidopsis_Evaluacion, by = c("sacc")) %>% dplyr::select(,-(9:12)) 
#resultados_arabidopsis_Evaluacion <- resultados_arabidopsis_Evaluacion %>% subset(go_id != "")
#####################################

#CUCUMIS-VITIS VINIFERA genes no anotados
df_CMVitis_all <- read.table("uva_all.txt",sep='\t', header = FALSE, col.names = nombres_columna)
filtrado_CMVitis_all <- df_CMVitis_all %>% filter(evalue <= MAX_EV)%>% filter(pidentity>=MIN_PID)%>% arrange(evalue)
filtrado_CMVitis_all <- filtrado_CMVitis_all %>% group_by(qacc) %>% distinct(qacc, .keep_all = TRUE) #nos quedamos con el resultado del gen con menor evalue

resultados_vitis <- data.frame("ensembl_gene_id"=c(""), 'ensembl_transcript_id'=c(""), "go_id"=c(""),"description"=c(""),"peptide"=c(""))

for (i in 1:nrow(filtrado_CMVitis_all)) {
  genCucumis <- filtrado_CMVitis_all[i, "qacc"]
  genVitis <- filtrado_CMVitis_all[i, "sacc"]
  
  matches_vitis <- list(genCucumis = list(genVitis))
  for (genCucumis in matches_vitis){
    genVitis <- matches_vitis[[1]][[1]]
    resultadosgenes_vitis <- biomaRt::getBM(attributes = c("ensembl_gene_id",'ensembl_transcript_id',"go_id","description","peptide"),mart= ensembl_vitis, filters = 'ensembl_transcript_id',values = unlist(genVitis))
    resultados_vitis <- full_join(resultados_vitis,resultadosgenes_vitis, by = c("ensembl_gene_id",'ensembl_transcript_id',"go_id","description","peptide"))
  }
}
resultados_vitis <- resultados_vitis[-c(1),]
colnames(resultados_vitis)[2]  <- "sacc"
resultados_vitis <- resultados_vitis %>% filter(go_id != "NA")%>% left_join(filtrado_CMVitis_all, by = c("sacc")) %>% dplyr::select(,-(9:12))
#####################################
#Para evaluación
df_CMVitis_Evaluacion <- read.table("uva_ganotados.txt",sep='\t', header = FALSE, col.names = nombres_columna)
filtrado_Vitis_Evaluacion <- df_CMVitis_Evaluacion %>% arrange(evalue) %>% filter(evalue <= MAX_EV) %>% distinct(qacc, .keep_all = TRUE) %>% filter(pidentity>=ID50) %>% filter(qcovs>=COV50)

resultados_vitis_Evaluacion <- data.frame("ensembl_gene_id"=c(""), 'ensembl_transcript_id'=c(""), "go_id"=c(""),"description"=c(""),"peptide"=c(""))

for (i in 1:nrow(filtrado_Vitis_Evaluacion)) {
  genCucumis <- filtrado_Vitis_Evaluacion[i, "qacc"]
  genVitis <- filtrado_Vitis_Evaluacion[i, "sacc"]
  
  matches_vitis <- list(genCucumis = list(genVitis))
  for (genCucumis in matches_vitis){
    genVitis <- matches_vitis[[1]][[1]]
    resultadosgenes_vitis <- biomaRt::getBM(attributes = c("ensembl_gene_id",'ensembl_transcript_id',"go_id","description","peptide"),mart= ensembl_vitis, filters = 'ensembl_transcript_id',values = unlist(genVitis))
    resultados_vitis_Evaluacion <- full_join(resultados_vitis_Evaluacion,resultadosgenes_vitis, by = c("ensembl_gene_id",'ensembl_transcript_id',"go_id","description","peptide"))
  }
}
resultados_vitis_Evaluacion <- resultados_vitis_Evaluacion[-c(1),]
colnames(resultados_vitis_Evaluacion)[2]  <- "sacc"
resultados_vitis_Evaluacion <- resultados_vitis_Evaluacion %>% filter(go_id != "NA")%>% left_join(filtrado_Vitis_Evaluacion, by = c("sacc")) %>% dplyr::select(,-(9:12))

#####################################

#Pasamos dataframes a lista de lista para la selección
colnames(gnoanotados_cucumis)[1]  <- "gen"
genesNoAnotados <- unique(as.list(gnoanotados_cucumis$gen)) #lista con los genes no anotados

#Arabidopsis thaliana - obtener la lista de listas con los resultados
anotacion_arabidopsis <- list()
for (i in 1:nrow(resultados_arabidopsis)) {
  gen_cucumis <- resultados_arabidopsis[i, "qacc"]
  go_id_ta <- resultados_arabidopsis[i, "go_id"]

  if (gen_cucumis %in% anotacion_arabidopsis){
    anotacion_arabidopsis[[gen_cucumis]] <- c(anotacion_arabidopsis[[gen_cucumis]],str_split(go_id_ta,";"))
    anotacion_arabidopsis <- anotacion_arabidopsis[[gen_cucumis]]
  } else {
    anotacion_arabidopsis <- c(anotacion_arabidopsis, list(gen_cucumis = str_split(go_id_ta,";")))
  }
}
names(anotacion_arabidopsis) <- resultados_arabidopsis$qacc
#####################################
#Para evaluación
anotacion_arabidopsis_Evaluacion <- list()

for (i in 1:nrow(resultados_arabidopsis_Evaluacion)) {
  gen_cucumis <- resultados_arabidopsis_Evaluacion[i, "qacc"]
  go_id_ta <- resultados_arabidopsis_Evaluacion[i, "go_id"]

  if (gen_cucumis %in% anotacion_arabidopsis_Evaluacion){
    anotacion_arabidopsis_Evaluacion[[gen_cucumis]] <- c(anotacion_arabidopsis_Evaluacion[[gen_cucumis]],str_split(go_id_ta,";"))
    anotacion_arabidopsis_Evaluacion <- anotacion_arabidopsis_Evaluacion[[gen_cucumis]]
  } else {
    anotacion_arabidopsis_Evaluacion <- c(anotacion_arabidopsis_Evaluacion, list(gen_cucumis = str_split(go_id_ta,";")))
  }
}
names(anotacion_arabidopsis_Evaluacion) <- resultados_arabidopsis_Evaluacion$qacc
#####################################
#Vitis vinifera - obtener la lista de listas con los resultados
anotacion_vitis <- list()
for (i in 1:nrow(resultados_vitis)) {
  gen_cucumis <- resultados_vitis[i, "qacc"]
  go_id_vitis <- resultados_vitis[i, "go_id"]

  if (gen_cucumis %in% resultados_vitis){
    anotacion_vitis[[gen_cucumis]] <- c(anotacion_vitis[[gen_cucumis]],str_split(go_id_vitis,";"))
    anotacion_vitis<- anotacion_vitis[[gen_cucumis]]
  } else {
    anotacion_vitis <- c(anotacion_vitis, list(gen_cucumis = str_split(go_id_vitis,";")))
  }
}
names(anotacion_vitis) <- resultados_vitis$qacc
#####################################
#Para evaluación
anotacion_vitis_Evaluacion <- list()

for (i in 1:nrow(resultados_vitis_Evaluacion)) {
  gen_cucumis <- resultados_vitis_Evaluacion[i, "qacc"]
  go_id_vitis <- resultados_vitis_Evaluacion[i, "go_id"]

  if (gen_cucumis %in% resultados_vitis_Evaluacion){
    anotacion_vitis_Evaluacion[[gen_cucumis]] <- c(anotacion_vitis_Evaluacion[[gen_cucumis]],str_split(go_id_vitis,";"))
    anotacion_vitis_Evaluacion<- anotacion_vitis_Evaluacion[[gen_cucumis]]
  } else {
    anotacion_vitis_Evaluacion <- c(anotacion_vitis_Evaluacion, list(gen_cucumis = str_split(go_id_vitis,";")))
  }
}

names(anotacion_vitis_Evaluacion) <- resultados_vitis_Evaluacion$qacc
#####################################
#Selección de anotaciones
una_especie <- list()
dos_especies <- list()
juntos <- list()

for (gen in genesNoAnotados){
  aTair <- list()
  aVitis <- list()
  if (gen %in% names(anotacion_arabidopsis)){
    aTair <- c(aTair,anotacion_arabidopsis[gen],use.names = TRUE)
    
  }
  if (gen %in% names(anotacion_vitis)){
    aVitis <-  c(aVitis, anotacion_vitis[gen],use.names = TRUE)
  }
  
  aP <- list()
  aP <- c(aTair,aVitis,use.names=TRUE)
  aP <- unlist(aP,recursive = FALSE)
  #print(aP)
  #names(aP_unique) <- names(aP)[match(aP_unique,aP)]
  #print(aP_unique)
  
  if (length(aP)>0){
    #Si solo encuentra anotaciones en un genoma
    if (length(aTair)==0 || length(aVitis)==0){
      una_especie <- c(una_especie,aP,use.names = TRUE)
    } else {
      for (a in names(aP)){
        if (a %in% names(aTair) && a %in% names(aVitis) && !(a %in% juntos)){
          comp1 <- aTair[a]
          comp2 <- aVitis[a]
          valores_coincidentes <- list(intersect(unlist(comp1), unlist(comp2)))
          names(valores_coincidentes) <- names(aVitis)
          dos_especies <- c(dos_especies,valores_coincidentes,use.names=TRUE)
          
          
          juntos <- c(juntos,names(comp1))
          juntos <- c(juntos,names(comp2))
          juntos <- unique(juntos)
          juntos <- c(juntos,juntos)
          
        }
      }
    }
  }
}

union <- c(una_especie,dos_especies,use.names = TRUE)
anotaciones_encontradasAM3 <- Filter(function(x) length(x)!= 0, lapply(union, function(x) x))
anotaciones_encontradasCH2 <- Filter(function(x) length(x)!= 0, lapply(union, function(x) x))
anotaciones_encontradasPSA2 <- Filter(function(x) length(x)!= 0, lapply(union, function(x) x))
anotaciones_encontradasPSU2 <- Filter(function(x) length(x)!= 0, lapply(union, function(x) x))
anotaciones_encontradasROCH2 <- Filter(function(x) length(x)!= 0, lapply(union, function(x) x))
anotaciones_encontradasTB2 <- Filter(function(x) length(x)!= 0, lapply(union, function(x) x))
anotaciones_encontradasTN2 <- Filter(function(x) length(x)!= 0, lapply(union, function(x) x))

#evaluación de anotación
colnames(ganotados_cucumis)[1] <-"gen" #Cambio nombre columna

genesAnotados <- list() #Lista con los genes anotados

for (i in 1:nrow(ganotados_cucumis)) {
  gen_cucumis <- ganotados_cucumis[i, "gen"]
  go_id_gnotados <- ganotados_cucumis[i, "go_id"]
  if (gen_cucumis %in% ganotados_cucumis){
    genesAnotados[[gen_cucumis]] <- c(genesAnotados[[gen_cucumis]],str_split(go_id_gnotados,";"))
    genesAnotados <- genesAnotados[[gen_cucumis]]
  } else {
    genesAnotados <- c(genesAnotados, list(gen_cucumis = str_split(go_id_gnotados,";")))
  }
}

names(genesAnotados) <- ganotados_cucumis$gen

anot_nopropuestas <- 0 #anotaciones en genes anotados que no se encuentran entre anotaciones propuestas
false_negative <- 0 #anotaciones que se encuentran en los genes anotados pero no en la predicción (FN)
false_positive <- 0 #anotaciones que se encuentran en predicción pero no en genes anotados (FP)
true_positive <- 0 #anotaciones que se encuentran en genes anotados y predicción (TP)

l <- list()
lb <- list()
juntos <- list()
for (gen in names(genesAnotados)) {
  aTair_e <- list()
  aVitis_e <- list()
  aSD <- list()
  
  aSD <- c(aSD,genesAnotados[gen],use.names = TRUE)
  
  if (gen %in% names(anotacion_arabidopsis_Evaluacion)){
    aTair_e <- c(aTair_e,anotacion_arabidopsis_Evaluacion[gen],use.names = TRUE)
  }
  
  if (gen %in% names(anotacion_vitis_Evaluacion)){
    aVitis_e <- c(aVitis_e,anotacion_vitis_Evaluacion[gen],use.names = TRUE)
  }
  
  aP <- list()
  aP <- c(aTair_e,aVitis_e,use.names=TRUE)
  aP <- unlist(aP,recursive = FALSE)
  #Comprobamos que anotaciones no coinciden con las obtenidas mediante blast
  for (i in aSD){
    if(!is.na(match(i, unlist(aVitis_e[gen]))) || (!is.na(match(i, unlist(aTair_e[gen]))))){
      anot_nopropuestas <- anot_nopropuestas + 1
    }
  }
  
  if (length(aP)>0){
    #Si solo encuentra anotaciones en un genoma
    if (length(aTair_e)==0 || length(aVitis_e)==0){
      l <- c(l,aP,use.names = TRUE)
    } else {
      for (a in names(aP)){
        if (a %in% names(aTair_e) && a %in% names(aVitis_e) && !(a %in% juntos)){
          comp1 <- aTair_e[a]
          comp2 <- aVitis_e[a]
          valores_coincidentes <- list(intersect(unlist(comp1), unlist(comp2))) #obtener la anotación que se repite en ambos organismos
          names(valores_coincidentes) <- names(aVitis_e)
          lb <- c(lb,valores_coincidentes,use.names=TRUE)
          
          juntos <- c(juntos,names(comp1)) #generamos un vector para que solo se haga el proceso 1 vez en un gen
          juntos <- c(juntos,names(comp2))
          juntos <- unique(juntos)
          juntos <- c(juntos,juntos)
        }
      }
    } 
  }
  union_e <- c(l,lb,use.names = TRUE)
  
  for (a in aSD){
    for (i in a){
      for(j in i){
        cat("Buscando: ",j)
        
        if (!is.na(match(j, unlist(union_e[gen])))) {
          true_positive <- true_positive +1
        } else {
          if (is.na(match(j, unlist(union_e[gen])))) {
            false_negative <- false_negative + 1
          }
        }
      }
    }
  }
  
  for (u in union_e[gen]){
    for (a in u){
      if (is.na(match(a, unlist(aSD)))){
        false_positive <- false_positive + 1
      }
    }
  }
}

cat("no encontradas: ", anot_nopropuestas)
cat("FN: ", false_negative)
cat("FP: ", false_positive)
cat("TP: ", true_positive)
cat("número de genes para los que se ha encontrado anotación: ", length(union))

precision <- ((true_positive/(true_positive+false_positive))*100)
print(precision)
exhaustividad <- (true_positive/(true_positive+false_negative)*100)
print(exhaustividad)
F1_score <- (((2*true_positive)/(2*true_positive+false_positive+false_negative))*100)
print(F1_score)
FDR <- (false_positive/(false_positive+true_positive))*100
print(FDR)

#################################   ANEXO   #########################

#gráficas para determinar rango de coverage
MAX_EV <- 0.000001
#para 75 genes
df_CMTomato <- df_CMTomato %>% filter(evalue <= MAX_EV) %>% arrange(evalue) %>% distinct(qacc, .keep_all = TRUE)
df_CMAmborella <- df_CMAmborella %>% filter(evalue <= MAX_EV) %>% arrange(evalue) %>% distinct(qacc, .keep_all = TRUE)
df_CMArabidopsis <- df_CMArabidopsis %>% filter(evalue <= MAX_EV) %>% arrange(evalue) %>% distinct(qacc, .keep_all = TRUE)
df_CMVitis <- df_CMVitis %>% filter(evalue <= MAX_EV) %>% arrange(evalue) %>% distinct(qacc, .keep_all = TRUE)

df_CMTomato$especie <- "Tomato"
df_CMAmborella$especie <- "Amborella"
df_CMVitis$especie <- "Vitis"
df_CMArabidopsis$especie <- "Arabidopsis"

join_resultados <-  merge(df_CMTomato,df_CMAmborella, all=TRUE)
join_resultados <-  merge(join_resultados,df_CMArabidopsis, all=TRUE)
join_resultados <-  merge(join_resultados,df_CMVitis, all=TRUE)

gf_coverage <- join_resultados %>% ggplot(aes(x=especie,y=qcovs)) + geom_violin(colour="black", fill="yellow")
gf_id <-join_resultados %>% ggplot(aes(x=especie,y=pidentity)) + geom_violin(colour="black", fill="blue")

#Dentro de un mismo rango de %id ¿son los mismos genes?
gTomato_rango <- df_CMTomato %>% filter(pidentity>=50 & pidentity<=60) %>% distinct(qacc, .keep_all = TRUE)
gArabidopsis_rango <- df_CMArabidopsis %>% filter(pidentity>=50 & pidentity<=60) %>% distinct(qacc, .keep_all = TRUE)
gAmborella_rango <- df_CMAmborella %>% filter(pidentity>=50 & pidentity<=60) %>% distinct(qacc, .keep_all = TRUE)
gVitis_rango <- df_CMVitis %>% filter(pidentity>=50 & pidentity<=60) %>% distinct(qacc, .keep_all = TRUE)

coincidencias1 <- gTomato_rango[ which(gTomato_rango$qacc %in% gArabidopsis_rango$qacc), ]
coincidencias2 <- gTomato_rango[ which(gTomato_rango$qacc %in% gAmborella_rango$qacc), ]
coincidencias3 <- gTomato_rango[ which(gTomato_rango$qacc %in% gVitis_rango$qacc), ]