library(tidyverse)
library(biomaRt)
library(openxlsx)

setwd("C:/Users/Andrea/Desktop/CUARTO/TFG/metilacion")

# Se juntan todos los data frame en una lista para facilitar la manipulación de los datos 
files <- Sys.glob("*.tsv") #Vector con el NOMBRE de los archivos de cada variedad

met_vector <- lapply(files, read.delim) #Lista con todos los df juntos
names_vector <- c()

for (file in files) { #Creación de un vector al que se le cambia el nombre del archivo para hacerlo más accesible
  name <- str_replace_all(file, "methylation_info_", "")
  name <- str_replace_all(name, ".tsv", "")
  names_vector[(length(names_vector) + 1)] = name #como un index
}

names(met_vector) <- names_vector #Asignación de nombres a cada df

for (i in 1:length(met_vector)) { #Adición de una nueva columna en la que se indique la variedad a la que pertenece
  df <- met_vector[[i]]
  df$variedad = names_vector[i] 
  met_vector[[i]] <- df
}

for (i in 1:length(met_vector)) {
  df <- met_vector[[i]]
  df <- df %>% dplyr::select(!(15:20)) #se eliminan columnas vacías
  met_vector[[i]] <- df
}

#Metilaciones totales de cada especie
for (i in 1:length(met_vector)) {
  df <- met_vector[[i]]
  print(count(df))
  met_vector[[i]] <- df
}
#Metilaciones totales por cromosoma
for (i in 1:length(met_vector)) {
  df <- met_vector[[i]]
  df <- df %>% filter(chr=="chr10") 
  met_vector[[i]] <- df
}

#Joins de los distintos data frames y coincidencias de metilación entre especies
join_am3 <- met_vector$`AM3-6` %>% left_join(met_vector$CH2,by= c("chr", "unique.id")) %>% left_join(met_vector$PSA2,by= c("chr", "unique.id")) %>% left_join(met_vector$PSU2,by= c("chr", "unique.id")) %>% left_join(met_vector$ROCH2,by= c("chr", "unique.id")) %>% left_join(met_vector$TB2,by= c("chr", "unique.id")) %>% left_join(met_vector$TN2,by= c("chr", "unique.id"))
num_igual_am3ch2 <- join_am3 %>% filter(join_am3[,28] != "NA")
num_igual_am3psa2 <- join_am3 %>% filter(join_am3[,41] != "NA")
num_igual_am3psu2 <- join_am3 %>% filter(join_am3[,54] != "NA")
num_igual_am3roch2 <- join_am3 %>% filter(join_am3[,67] != "NA")
num_igual_am3tb2 <- join_am3 %>% filter(join_am3[,80] != "NA")
num_igual_am3tn2 <- join_am3 %>% filter(join_am3[,93] != "NA")

join_ch2 <- met_vector$CH2 %>% left_join(met_vector$`AM3-6`,by= c("chr", "unique.id")) %>% left_join(met_vector$PSA2,by= c("chr", "unique.id")) %>% left_join(met_vector$PSU2,by= c("chr", "unique.id")) %>% left_join(met_vector$ROCH2,by= c("chr", "unique.id")) %>% left_join(met_vector$TB2,by= c("chr", "unique.id")) %>% left_join(met_vector$TN2,by= c("chr", "unique.id"))
num_igual_ch2psa2 <- join_ch2 %>% filter(join_ch2[,41] != "NA")
num_igual_ch2psu2 <- join_ch2 %>% filter(join_ch2[,54] != "NA")
num_igual_ch2roch2 <- join_ch2 %>% filter(join_ch2[,67] != "NA")
num_igual_ch2tb2 <- join_ch2 %>% filter(join_ch2[,80] != "NA")
num_igual_ch2tn2 <- join_ch2 %>% filter(join_ch2[,93] != "NA")

join_psa2 <- met_vector$PSA2 %>% left_join(met_vector$CH2,by= c("chr", "unique.id")) %>% left_join(met_vector$`AM3-6`,by= c("chr", "unique.id")) %>% left_join(met_vector$PSU2,by= c("chr", "unique.id")) %>% left_join(met_vector$ROCH2,by= c("chr", "unique.id")) %>% left_join(met_vector$TB2,by= c("chr", "unique.id")) %>% left_join(met_vector$TN2,by= c("chr", "unique.id"))
num_igual_psa2psu2 <- join_psa2 %>% filter(join_psa2[,54] != "NA")
num_igual_psa2roch2 <- join_psa2 %>% filter(join_psa2[,67] != "NA")
num_igual_psa2tb2 <- join_psa2 %>% filter(join_psa2[,80] != "NA")
num_igual_psa2tn2 <- join_psa2 %>% filter(join_psa2[,93] != "NA")

join_psu2 <- met_vector$PSU2 %>% left_join(met_vector$CH2,by= c("chr", "unique.id")) %>% left_join(met_vector$`AM3-6`,by= c("chr", "unique.id")) %>% left_join(met_vector$PSA2,by= c("chr", "unique.id")) %>% left_join(met_vector$ROCH2,by= c("chr", "unique.id")) %>% left_join(met_vector$TB2,by= c("chr", "unique.id")) %>% left_join(met_vector$TN2,by= c("chr", "unique.id"))
num_igual_psu2roch2 <- join_psu2 %>% filter(join_psu2[,67] != "NA")
num_igual_psu2tb2 <- join_psu2 %>% filter(join_psu2[,80] != "NA")
num_igual_psu2tn2 <- join_psu2 %>% filter(join_psu2[,93] != "NA")

join_roch2 <- met_vector$ROCH2 %>% left_join(met_vector$CH2,by= c("chr", "unique.id")) %>% left_join(met_vector$`AM3-6`,by= c("chr", "unique.id")) %>% left_join(met_vector$PSA2,by= c("chr", "unique.id")) %>% left_join(met_vector$PSU2,by= c("chr", "unique.id")) %>% left_join(met_vector$TB2,by= c("chr", "unique.id")) %>% left_join(met_vector$TN2,by= c("chr", "unique.id"))
num_igual_roch2tb2 <- join_roch2 %>% filter(join_roch2[,80] != "NA")
num_igual_roch2tn2 <- join_roch2 %>% filter(join_roch2[,93] != "NA")

join_tb2 <- met_vector$TB2 %>% left_join(met_vector$CH2,by= c("chr", "unique.id")) %>% left_join(met_vector$`AM3-6`,by= c("chr", "unique.id")) %>% left_join(met_vector$PSA2,by= c("chr", "unique.id")) %>% left_join(met_vector$PSU2,by= c("chr", "unique.id")) %>% left_join(met_vector$ROCH2,by= c("chr", "unique.id")) %>% left_join(met_vector$TN2,by= c("chr", "unique.id"))
num_igual_tb2tn2 <- join_tb2 %>% filter(join_tb2[,93] != "NA")

join_tn2 <- met_vector$TN2 %>% left_join(met_vector$CH2,by= c("chr", "unique.id")) %>% left_join(met_vector$`AM3-6`,by= c("chr", "unique.id")) %>% left_join(met_vector$PSA2,by= c("chr", "unique.id")) %>% left_join(met_vector$PSU2,by= c("chr", "unique.id")) %>% left_join(met_vector$ROCH2,by= c("chr", "unique.id")) %>% left_join(met_vector$TB2,by= c("chr", "unique.id"))

#Metilaciones generales pertenecientes al pangenoma
delete.na <- function(df, n=0) {
  df[rowSums(is.na(df)) <= n,]
}
pangenoma_metilacion <- delete.na(join_am3)

p_chr00 <- pangenoma_metilacion %>% filter(chr=="chr00")%>% filter(Annotation.x !="Intergenic")
p_chr01 <- pangenoma_metilacion %>% filter(chr=="chr01")%>% filter(Annotation.x !="Intergenic")
p_chr02 <- pangenoma_metilacion %>% filter(chr=="chr02")%>% filter(Annotation.x !="Intergenic")
p_chr03 <- pangenoma_metilacion %>% filter(chr=="chr03")%>% filter(Annotation.x !="Intergenic")
p_chr04 <- pangenoma_metilacion %>% filter(chr=="chr04")%>% filter(Annotation.x !="Intergenic")
p_chr05 <- pangenoma_metilacion %>% filter(chr=="chr05")%>% filter(Annotation.x !="Intergenic")
p_chr06 <- pangenoma_metilacion %>% filter(chr=="chr06")%>% filter(Annotation.x !="Intergenic")
p_chr07 <- pangenoma_metilacion %>% filter(chr=="chr07")%>% filter(Annotation.x !="Intergenic")
p_chr08 <- pangenoma_metilacion %>% filter(chr=="chr08")%>% filter(Annotation.x !="Intergenic")
p_chr09 <- pangenoma_metilacion %>% filter(chr=="chr09")%>% filter(Annotation.x !="Intergenic")
p_chr10 <- pangenoma_metilacion %>% filter(chr=="chr10")%>% filter(Annotation.x !="Intergenic")
p_chr11 <- pangenoma_metilacion %>% filter(chr=="chr11")%>% filter(Annotation.x !="Intergenic")
p_chr12 <- pangenoma_metilacion %>% filter(chr=="chr12")%>% filter(Annotation.x !="Intergenic")
 
#Metilaciones únicas de cada variante
metilacionesunicas_am3 <- join_am3 %>% filter(is.na(variedad.y)) %>% filter(is.na(variedad.x.x)) %>% filter(is.na(variedad.y.y)) %>%filter(is.na(variedad.x.x.x)) %>% filter(is.na(variedad.y.y.y)) %>% filter(is.na(variedad)) %>% filter(Annotation.x !="Intergenic")
metilacionesunicas_ch2 <- join_ch2 %>% filter(is.na(variedad.y)) %>% filter(is.na(variedad.x.x)) %>% filter(is.na(variedad.y.y)) %>%filter(is.na(variedad.x.x.x)) %>% filter(is.na(variedad.y.y.y)) %>% filter(is.na(variedad)) %>% filter(Annotation.x !="Intergenic")
metilacionesunicas_psa2 <- join_psa2 %>% filter(is.na(variedad.y)) %>% filter(is.na(variedad.x.x)) %>% filter(is.na(variedad.y.y)) %>%filter(is.na(variedad.x.x.x)) %>% filter(is.na(variedad.y.y.y)) %>% filter(is.na(variedad)) %>% filter(Annotation.x !="Intergenic")
metilacionesunicas_psu2 <- join_psu2 %>% filter(is.na(variedad.y)) %>% filter(is.na(variedad.x.x)) %>% filter(is.na(variedad.y.y)) %>%filter(is.na(variedad.x.x.x)) %>% filter(is.na(variedad.y.y.y)) %>% filter(is.na(variedad)) %>% filter(Annotation.x !="Intergenic")
metilacionesunicas_roch2 <- join_roch2 %>% filter(is.na(variedad.y)) %>% filter(is.na(variedad.x.x)) %>% filter(is.na(variedad.y.y)) %>%filter(is.na(variedad.x.x.x)) %>% filter(is.na(variedad.y.y.y)) %>% filter(is.na(variedad)) %>% filter(Annotation.x !="Intergenic")
metilacionesunicas_tb2 <- join_tb2 %>% filter(is.na(variedad.y)) %>% filter(is.na(variedad.x.x)) %>% filter(is.na(variedad.y.y)) %>%filter(is.na(variedad.x.x.x)) %>% filter(is.na(variedad.y.y.y)) %>% filter(is.na(variedad)) %>% filter(Annotation.x !="Intergenic")
metilacionesunicas_tn2 <- join_tn2 %>% filter(is.na(variedad.y)) %>% filter(is.na(variedad.x.x)) %>% filter(is.na(variedad.y.y)) %>%filter(is.na(variedad.x.x.x)) %>% filter(is.na(variedad.y.y.y)) %>% filter(is.na(variedad)) %>% filter(Annotation.x !="Intergenic")

#Obtención de nombres y anotación de genes con biomart

ensembl_plants <- useEnsemblGenomes(biomart = "plants_mart")
searchDatasets(ensembl_plants, pattern = "Cucumis melo")
ensembl_cucumis <- useEnsemblGenomes(biomart = "plants_mart", 
                                     dataset = "cmelo_eg_gene")

#Anotación cucumis
anotacion_cucumis <- biomaRt::getBM(attributes = c("ensembl_gene_id","go_id","description","peptide"),mart= ensembl_cucumis)
colnames(anotacion_cucumis)[3]  <- "Nearest.Unigene.x" 

#Genes totales metilados anotados / no anotados de cada variedad
unique.genes <- function(df){ #conocer el número de genes metilados en cada variedad
  df <- df %>% distinct(Nearest.Unigene.x, .keep_all = TRUE) %>% dplyr::select(!(16:93)) 
}

metilacionesunicas_am3 <- unique.genes(metilacionesunicas_am3)
metilacionesunicas_ch2 <- unique.genes(metilacionesunicas_ch2)
metilacionesunicas_psa2 <- unique.genes(metilacionesunicas_psa2)
metilacionesunicas_psu2 <- unique.genes(metilacionesunicas_psu2)
metilacionesunicas_roch2 <- unique.genes(metilacionesunicas_roch2)
metilacionesunicas_tb2 <- unique.genes(metilacionesunicas_tb2)
metilacionesunicas_tn2 <- unique.genes(metilacionesunicas_tn2)

############################### obtener aquellos genes que carecen de anotación para el script "tBlast_Annotation.R" ###############################
annotation_search <- function(df){ 
  ganotados_cucumis <- merge(df,anotacion_cucumis,by= c("Nearest.Unigene.x")) 
  ganotados_cucumis <- ganotados_cucumis %>% dplyr::select("Nearest.Unigene.x","chr","start.x","stop.x","Annotation.x","description","go_id","peptide")
  
  gnoanotados <- ganotados_cucumis %>% filter(go_id=="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)
}

gnoanotados_cucumis <- annotation_search(metilacionesunicas_am3)
gnoanotados_cucumis <- annotation_search(metilacionesunicas_ch2)
gnoanotados_cucumis <- annotation_search(metilacionesunicas_psa2)
gnoanotados_cucumis <- annotation_search(metilacionesunicas_psu2)
gnoanotados_cucumis <- annotation_search(metilacionesunicas_roch2)
gnoanotados_cucumis <- annotation_search(metilacionesunicas_tb2)
gnoanotados_cucumis <- annotation_search(metilacionesunicas_tn2)
 ###################################################################################################################################################

#AM3-6
ganotados_am3 <- merge(metilacionesunicas_am3,anotacion_cucumis, by= c("Nearest.Unigene.x"))
ganotados_am3 <- ganotados_am3 %>% dplyr::select("Nearest.Unigene.x","chr","start.x","stop.x","Annotation.x","description","go_id","peptide")
  
gnoanotados_am3 <- ganotados_am3 %>% filter(go_id=="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)
ganotados_am3 <- ganotados_am3 %>% filter(go_id!="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)

  #CH2
ganotados_ch2 <- merge(metilacionesunicas_ch2,anotacion_cucumis, by= c("Nearest.Unigene.x"))
ganotados_ch2 <- ganotados_ch2 %>% dplyr::select("Nearest.Unigene.x","chr","start.x","stop.x","Annotation.x","description","go_id","peptide")

gnoanotados_ch2 <- ganotados_ch2 %>% filter(go_id=="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)
ganotados_ch2 <- ganotados_ch2 %>% filter(go_id!="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)

  #PSA2
ganotados_psa2 <- merge(metilacionesunicas_psa2,anotacion_cucumis, by= c("Nearest.Unigene.x"))
ganotados_psa2 <- ganotados_psa2 %>% dplyr::select("Nearest.Unigene.x","chr","start.x","stop.x","Annotation.x","description","go_id","peptide")

gnoanotados_psa2 <- ganotados_psa2 %>% filter(go_id=="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)
ganotados_psa2 <- ganotados_psa2 %>% filter(go_id!="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)

  #PSU2
ganotados_psu2 <- merge(metilacionesunicas_psu2,anotacion_cucumis, by= c("Nearest.Unigene.x"))
ganotados_psu2 <- ganotados_psu2 %>% dplyr::select("Nearest.Unigene.x","chr","start.x","stop.x","Annotation.x","description","go_id","peptide")

gnoanotados_psu2 <- ganotados_psu2 %>% filter(go_id=="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)
ganotados_psu2 <- ganotados_psu2 %>% filter(go_id!="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)

  #ROCH2
ganotados_roch2 <- merge(metilacionesunicas_roch2,anotacion_cucumis, by= c("Nearest.Unigene.x"))
ganotados_roch2 <- ganotados_roch2 %>% dplyr::select("Nearest.Unigene.x","chr","start.x","stop.x","Annotation.x","description","go_id","peptide")

gnoanotados_roch2 <- ganotados_roch2 %>% filter(go_id=="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)
ganotados_roch2 <- ganotados_roch2 %>% filter(go_id!="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)

  #TB2
ganotados_tb2 <- merge(metilacionesunicas_tb2,anotacion_cucumis, by= c("Nearest.Unigene.x"))
ganotados_tb2 <- ganotados_tb2 %>% dplyr::select("Nearest.Unigene.x","chr","start.x","stop.x","Annotation.x","description","go_id","peptide")

gnoanotados_tb2 <- ganotados_tb2 %>% filter(go_id=="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)
ganotados_tb2 <- ganotados_tb2 %>% filter(go_id!="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)

  #TN2
ganotados_tn2 <- merge(metilacionesunicas_tn2,anotacion_cucumis, by= c("Nearest.Unigene.x"))
ganotados_tn2 <- ganotados_tn2 %>% dplyr::select("Nearest.Unigene.x","chr","start.x","stop.x","Annotation.x","description","go_id","peptide")

gnoanotados_tn2 <- ganotados_tn2 %>% filter(go_id=="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)
ganotados_tn2 <- ganotados_tn2 %>% filter(go_id!="") %>% distinct(Nearest.Unigene.x, .keep_all = TRUE)



###########################         ANEXO           ###################################

#Obtener los genes con metilaciones bases/pangenoma metilación
lista_cromosomas <- list(p_chr00,p_chr01,p_chr02,p_chr03,p_chr04,p_chr05,p_chr06,p_chr07,p_chr08,p_chr09,p_chr10,p_chr11,p_chr12)
genes <- c("start.x","stop.x","Annotation.x")
resultados_genes_metilacionbase <- lapply(lista_cromosomas, "[", genes)
print(resultados_genes_metilacionbase[[1]][3])
write.xlsx(resultados_genes_metilacionbase[[14]][3],"genes_ch13.xlsx")
#Obtener el ID de los genes metilados en distintas variedades en excel
for (i in 1:length(met_vector)) {
  df <- met_vector[[i]]
  df <- df %>% filter(Annotation !="Intergenic") #se eliminan las metilaciones intergénicas
  met_vector[[i]] <- df
}
genes <- c("chr","start","stop","Annotation")
resultados_genes <-  lapply(met_vector, "[", genes)

write.xlsx(resultados_genes$`AM3-6`,"AM3genes.xlsx")
write.xlsx(resultados_genes$CH2,"CH2genes.xlsx")
write.xlsx(resultados_genes$PSA2,"PSA2genes.xlsx")
write.xlsx(resultados_genes$PSU2,"PSU2genes.xlsx")
write.xlsx(resultados_genes$ROCH2,"ROCHgenes.xlsx")
write.xlsx(resultados_genes$TB2,"TB2genes.xlsx")
write.xlsx(resultados_genes$TN2,"TN2genes.xlsx")

#Genes metilaciones únicas de cada variante en excel
lista_variedades <- list(metilacionesunicas_am3,metilacionesunicas_ch2,metilacionesunicas_psa2,metilacionesunicas_psu2,metilacionesunicas_roch2,metilacionesunicas_tb2,metilacionesunicas_tn2)
genes <- c("chr","start.x","stop.x","Annotation.x")
resultados_genes <- lapply(lista_variedades, "[", genes)

write.xlsx(resultados_genes[[1]],"AM3genesunicos.xlsx")
write.xlsx(resultados_genes[[2]],"CH2genesunicos.xlsx")
write.xlsx(resultados_genes[[3]],"PSA2genesunicos.xlsx")
write.xlsx(resultados_genes[[4]],"PSU2genesunicos.xlsx")
write.xlsx(resultados_genes[[5]],"ROCHgenesunicos.xlsx")
write.xlsx(resultados_genes[[6]],"TB2genesunicos.xlsx")
write.xlsx(resultados_genes[[7]],"TN2genesunicos.xlsx")

