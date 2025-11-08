## Base utilities for building and handling gene trees ##
## This file contains helper functions used by the Simulated annealing and NSGA2 code:
## - loading aligned sequences and computing pairwise distances
## - building neighbor-joining trees for a set of genes
## - convenience functions to initialize populations and read gene sets
## files where this is used: nsga2.r, simAnnMO.r, main.r, mutation.r, crossover.r

if( !exists("path") ) path = "~/Documents/U2021/"
library("phangorn")
library("readxl")
source(paste0(path,"code_base/used_packages.r"))

### Global variables ###
# directory holding aligned sequence files (phylip format expected)
directory = paste0(path, "Secuencias_alineadas/")
# list of filenames in 'directory' (used by index-based selection)
genesIndex = list.files(directory)
# reference tree used to score candidate trees (RF distance)
modelTree = read.tree(paste0(path,"Instances/arbol_referencia_nuevo.tree"))
# Note: genes are used according to their position in this list. Gene names are taken from here.

### Functions ###

buildTree <- function(genesList){
  # Build a neighbor-joining tree from a list of gene indices.
  # genesList: vector of integer indices into `genesIndex`.
  write.table(3, file = "out.txt", sep = "\n",
              row.names = TRUE, col.names = NA)
  mat = NULL
  for(i in 1:length(genesList)){
    # read the aligned sequence file for this gene index
    archivo <- read.dna(file=paste(directory, genesIndex[genesList[i]], sep=""), format = "sequential")
    datos <- phyDat(data=archivo, type = "DNA")
    # accumulate Hamming distances across genes (summing distance matrices)
    mat <- if(is.null(mat)) dist.hamming(datos) else mat + dist.hamming(datos)
  }
  arbol <- nj(mat) # neighbor-joining tree from the aggregated distance matrix
  arbol$genes = genesList
  # score tree using normalized Robinson-Foulds distance against modelTree
  arbol$score = RF.dist(arbol, tree2 = modelTree, normalize = TRUE)
  arbol$ranking = 0
  # arbol$flag retained from older experiments; left unused here
  return(arbol)
}

getRandomGenes <- function(size){
  # Return `size` random indices in 1:length(genesIndex)
  return(sample(1:length(genesIndex), size))
}

buildTrees <- function(genes_list){
  # Build multiple trees given a list of gene index vectors.
  population = NULL
  
  for(genes in genes_list){
    population[[i]] = buildTree(genes)
  }
  return(population)
}

init <- function(pop_size, max_genes){
  # Initialize a population of `pop_size` trees. Each tree is built from
  # a random number (1..max_genes) of randomly selected genes.
  population = NULL
  for(i in 1:pop_size){
    population[[i]] = buildTree(getRandomGenes(sample(1:max_genes, 1)))
  }
  return(population)
}

setNewGeneIndex <- function(newIndexNames){
  # Update the global `genesIndex` with new filenames (append file suffix)
  for(i in 1:length(newIndexNames)){
    newIndexNames[i] = paste0(" ", newIndexNames[i], ".fasta.phylip")
  }
  assign("genesIndex", newIndexNames, envir = .GlobalEnv)
}

changeModelTree <- function(file){
  # Replace the global modelTree with a new tree loaded from `file`.
  assign("modelTree", read.tree(file), envir = .GlobalEnv)
  #print(paste("Changed reference tree with: ",file, sep=""))
}

buildTreeFromXSLX <- function(file_name, genes_set){
  # Read an Excel sheet that describes gene groups / medoids and convert
  # the selection into indices used by buildTree.
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

  # convert gene names into the same filename format used by genesIndex
  for(i in 1:length(ref_genes)){
    ref_genes[i] = paste0(" ", ref_genes[i], ".fasta.phylip")
  }

  # match returns positions (indices) of ref_genes inside genesIndex
  return(buildTree(match(ref_genes, genesIndex)))
}
