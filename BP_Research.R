#se estudia la relación de los genes anotados con procesos biológicos
lista_ganotados <- list(total_ganotados_am3,total_ganotados_ch2,total_ganotados_psa2,total_ganotados_psu2,total_ganotados_roch2,total_ganotados_tb2,total_ganotados_tn2)
lista_variedades <- list("AM3","CH2","PSA2","PSU2","ROCH2","TB2","TN2")
names(lista_ganotados) <- lista_variedades

#terpenoides
go_terpenoides <- list("GO:0006714","GO:0016114")

analysis.terpenoides<- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_terpenoides)))){
            cat(anotacion,": ","sintesis de terpenoides. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.terpenoides(lista_ganotados[i])
}  


#indol
go_indol <- list("GO:0034768","GO:0042343","GO:0035834","GO:0042433","GO:0042432","GO:0033984")

analysis.indol<- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_indol)))){
            cat(anotacion,": ","sintesis de indol ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.indol(lista_ganotados[i])
}  

#metilcinamato
go_cinnamate <- list("GO:0016710","GO:0106290","GO:0050344","GO:0050412","GO:0042431","GO:0043786")

analysis.cinnamate <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_cinnamate)))){
            cat(anotacion,": ","sintesis de cinamato. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.cinnamate(lista_ganotados[i])
}  


#sesquiterpenos
go_sesquiterpenos <- list("GO:0016106","GO:0016107","GO:0051763","GO:0051761")

analysis.sesquiterpenos <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_sesquiterpenos)))){
            cat(anotacion,": ","sintesis de sesquiterpenos. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.sesquiterpenos(lista_ganotados[i])
}  

#sesquiterpenoides
go_sesquiterpenoides <- list("GO:0016106","GO:0016107")

analysis.sesquiterpenoides <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_sesquiterpenoides)))){
            cat(anotacion,": ","sintesis de sesquiterpenoides. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.sesquiterpenoides(lista_ganotados[i])
}  

#alcohol
go_alcohol <- list("GO:0018456","GO:0043178","GO:0006066","GO:0046445")
analysis.alcohol <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_alcohol)))){
            cat(anotacion,": ","sintesis de alcohol. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.alcohol(lista_ganotados[i])
}  

#bencilisoquinolina
go_bisoquinoline <- list("GO:0009708")
analysis.bisoquinoline <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_bisoquinoline)))){
            cat(anotacion,": ","sintesis de bencilisoquinolina. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.bisoquinoline(lista_ganotados[i])
}  

#bencilalcohol
go_balcohol <- list("GO:0102703")
analysis.balcohol <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_balcohol)))){
            cat(anotacion,": ","sintesis de bencilalcohol. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.balcohol(lista_ganotados[i])
}  

#canfeno
go_canfeno <- list("GO:0102703","GO:0043693")

analysis.canfeno <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_canfeno)))){
            cat(anotacion,": ","sintesis de canfeno. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.canfeno(lista_ganotados[i])
}  

#monoterpenos
go_monoterpenos <- list("GO:0043692","GO:0043694")
analysis.monoterpenos <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_monoterpenos)))){
            cat(anotacion,": ","sintesis de monoterpenos. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.monoterpenos(lista_ganotados[i])
} 

#ocimeno
go_ocimeno <- list("GO:0042408")
analysis.ocimeno <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_ocimeno)))){
            cat(anotacion,": ","sintesis de ocimenos. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.ocimeno(lista_ganotados[i])
} 

#careno
go_careno <- list("GO:0102702")
analysis.careno <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_careno)))){
            cat(anotacion,": ","sintesis de careno. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.careno(lista_ganotados[i])
} 

#cariofileno
go_cariofileno <- list("GO:0080016")

analysis.cariofileno <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_cariofileno)))){
            cat(anotacion,": ","sintesis de cariofileno. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.cariofileno(lista_ganotados[i])
} 

#limoneno
go_limoneno <- list("GO:0034002")

analysis.limoneno <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_limoneno)))){
            cat(anotacion,": ","sintesis de limoneno. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.limoneno(lista_ganotados[i])
} 

#metil benzoato
go_mbenzoato <- list("GO:0080150")

analysis.mbenzoato <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_mbenzoato)))){
            cat(anotacion,": ","sintesis de metilbenzoato. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.mbenzoato(lista_ganotados[i])
} 

#selineno
go_selineno <- list("GO:0102906")

analysis.selineno <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_selineno)))){
            cat(anotacion,": ","sintesis de selineno. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.selineno(lista_ganotados[i])
} 

#humuleno
go_humuleno<- list("GO:0080017","GO:0080016")

analysis.humuleno <- function(ganotados){
  for (listas in ganotados){
    for (gen in listas){
      for (anotaciones in gen) {
        for (anotacion in anotaciones){
          if(!is.na(match(anotacion, unlist(go_humuleno)))){
            cat(anotacion,": ","sintesis de humuleno. ")
          }
        }
      }
    }
  }
}

for (i in 1:length(lista_ganotados)) { 
  print(names(lista_ganotados[i]))
  analysis.humuleno(lista_ganotados[i])
} 

