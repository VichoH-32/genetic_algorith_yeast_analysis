if( !exists("path") ) path = "~/Documents/U2021/"
library("phangorn")
library("readxl")
source(paste0(path,"code_base/used_packages.r"))
### Global variables ###
directory = paste0(path, "Secuencias_alineadas/")
genesIndex = list.files(directory)
modelTree = read.tree(paste0(path,"Instances/arbol_referencia_nuevo.tree"))
# Nota: los genes se trabajan segun la posición del gen esta lista. De aquí se saca el nombre del gen.

### Functions ###

buildTree <- function(genesList){
  #print(genesList)
  write.table(3, file = "out.txt", sep = "\n",
              row.names = TRUE, col.names = NA)
  mat = NULL
  for(i in 1:length(genesList)){
    #print(paste("ite: ", i, " of ",length(genesList) , ", gen: ", genesIndex[genesList[i]], sep=""))
    archivo <- read.dna(file=paste(directory, genesIndex[genesList[i]], sep=""), format = "sequential")
    datos <- phyDat(data=archivo, type = "DNA")
    mat <- if(is.null(mat)) dist.hamming(datos) else mat + dist.hamming(datos)
  }
  arbol <- nj(mat)
  arbol$genes = genesList
  arbol$score = RF.dist(arbol, tree2 = modelTree, normalize = TRUE)
  arbol$ranking = 0
  #arbol$flag = FALSE # no recuerdo su objetivo XD
  
  return(arbol)
}

getRandomGenes <- function(size){
  return(sample(1:length(genesIndex), size)) # get from 1 to size random positions from genesIndex list
}

buildTrees <- function(genes_list){
  population = NULL
  
  for(genes in genes_list){
    population[[i]] = buildTree(genes)
  }
  return(population)
}

init <- function(pop_size, max_genes){
  population = NULL
  
  for(i in 1:pop_size){
    population[[i]] = buildTree(getRandomGenes(sample(1:max_genes, 1)))
  }
  return(population)
}

setNewGeneIndex <- function(newIndexNames){
  for(i in 1:length(newIndexNames)){
    newIndexNames[i] = paste0(" ", newIndexNames[i], ".fasta.phylip")
  }
  assign("genesIndex", newIndexNames, envir = .GlobalEnv)
}

changeModelTree <- function(file){
  assign("modelTree", read.tree(file), envir = .GlobalEnv)
  #print(paste("Changed reference tree with: ",file, sep=""))
}

buildTreeFromXSLX <- function(file_name, genes_set){
  df = read_excel(paste0(path, file_name))
  
  if(genes_set == 1){
    ref_genes = df[df$"GRUPOS_CON_REF" == 1 & df$GENES != "REF",]$GENES
  }else if(genes_set == 2){
    ref_genes = df[df$"MEDOIDE_CON_REF" == "SI" & df$GENES != "REF",]$GENES
  }else if(genes_set == 3){
    ref_genes = df[df$"GRUPOS_SIN_REF" == 1 & df$GENES != "REF",]$GENES
  }else if(genes_set == 4){
    ref_genes = df[df$"MEDOIDE_SIN_REF" == "SI" & df$GENES != "REF",]$GENES
  }else{
    ref_genes = df[df$"ARTICULO_ANTIGUO" == "SI" & df$GENES != "REF",]$GENES
  }
  
  for(i in 1:length(ref_genes)){
    ref_genes[i] = paste0(" ", ref_genes[i], ".fasta.phylip")
  }
  
  return(buildTree(match(ref_genes, genesIndex)))
}
