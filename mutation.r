if( !exists("path") ) path = "~/Documents/U2021/"
library("phangorn")
source(paste0(path,"code_base/base.r"))

### Functions ###

mutateGenes <- function(prev_genes, percentage){
  if(length(prev_genes) > 1){
    genesReplacePositions = sample(1:length(prev_genes), floor( length(prev_genes) * percentage))
  }else{
    genesReplacePositions = sample(1:length(prev_genes), 1)
  }
  genesReplacementsIndex = getRandomGenes(length(genesReplacePositions))
  
  # cambia los genes
  for(i in 1:length(genesReplacePositions)){prev_genes[genesReplacePositions[i]] = genesReplacementsIndex[i]}
  return(prev_genes)
}

mutateLength <- function(prev_genes, max_length){
  sl = sample(1:length(prev_genes), 1)
  prev_genes = unique(append( sample(prev_genes, sl), getRandomGenes(sample(1:length(prev_genes), 1))))
  
  if(length(prev_genes) > max_length){return(prev_genes[1:max_length])}
  return(prev_genes)
}

mutateTree <- function(prev_genes, percentage, can_mutate_length=FALSE, max_length){
  mutatedGenes = prev_genes
  if(can_mutate_length){
    if(missing(max_length)){
      stop("if can mutate length is enabled max_length can't be NULL")
    }
    mutatedGenes = mutateLength(mutatedGenes, max_length)
  }
  
  mutatedGenes = mutateGenes(mutatedGenes, percentage)
  
  # reconstruye el arbol con los genes nuevos (sin considerar la bandera)
  return(buildTree(mutatedGenes))
}

mutation <- function(population, percentage, can_mutate_length=FALSE, max_length){
  if(percentage > 1 | percentage < 0){stop("Given percentage is not valid")}
  
  for(i in 1:length(population)){ # recorre toda la poblacion
    # muta el arbol i
    if(can_mutate_length){
      if(missing(max_length)){
        stop("if can mutate length is enabled max_length can't be NULL")
      }
      population[[i]] = mutateTree(population[[i]]$genes, percentage, TRUE, max_genes)
    }else{
      population[[i]] = mutateTree(population[[i]]$genes, percentage)
    }
  }
  return(population)
}

init2 <- function(initial_solution, pop_size, mutation_percentage, max_genes){ # only used when an initial solution is provided
  population = NULL
  
  population[[1]] = initial_solution
  write.table(2, file = "out.txt", sep = "\n",
              row.names = TRUE, col.names = NA)
  for(i in 2:pop_size){
    population[[i]] = mutateTree(population[[1]]$genes, mutation_percentage, TRUE, max_genes)
  }
  return(population)
}