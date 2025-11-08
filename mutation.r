## Mutation operators and helpers ##
## Functions to mutate gene lists (change genes and optionally change length)

if( !exists("path") ) path = "~/Documents/U2021/"
library("phangorn")
source(paste0(path,"code_base/base.r"))

### Functions ###

mutateGenes <- function(prev_genes, percentage){
  # Replace a fraction (percentage) of genes in prev_genes with random genes.
  if(length(prev_genes) > 1){
    genesReplacePositions = sample(1:length(prev_genes), floor( length(prev_genes) * percentage))
  }else{
    # if only one gene exists, choose position 1 to replace
    genesReplacePositions = sample(1:length(prev_genes), 1)
  }
  genesReplacementsIndex = getRandomGenes(length(genesReplacePositions))
  
  # Replace selected positions with random gene indices
  for(i in 1:length(genesReplacePositions)){prev_genes[genesReplacePositions[i]] = genesReplacementsIndex[i]}
  return(prev_genes)
}

mutateLength <- function(prev_genes, max_length){
  # Randomly alter the length of gene list by sampling a subsequence and
  # appending a few random genes. Ensures uniqueness of genes.
  sl = sample(1:length(prev_genes), 1)
  prev_genes = unique(append( sample(prev_genes, sl), getRandomGenes(sample(1:length(prev_genes), 1))))
  
  if(length(prev_genes) > max_length){return(prev_genes[1:max_length])}
  return(prev_genes)
}

mutateTree <- function(prev_genes, percentage, can_mutate_length=FALSE, max_length){
  # Build a mutated gene list and then rebuild the tree for it.
  mutatedGenes = prev_genes
  if(can_mutate_length){
    if(missing(max_length)){
      stop("if can mutate length is enabled max_length can't be NULL")
    }
    mutatedGenes = mutateLength(mutatedGenes, max_length)
  }
  mutatedGenes = mutateGenes(mutatedGenes, percentage)
  # Reconstruct tree for the mutated gene set
  return(buildTree(mutatedGenes))
}

mutation <- function(population, percentage, can_mutate_length=FALSE, max_length){
  # Apply mutation to all individuals in `population`.
  if(percentage > 1 | percentage < 0){stop("Given percentage is not valid")}
  
  for(i in 1:length(population)){ # iterate population mutating tree "i"
    if(can_mutate_length){
      if(missing(max_length)){
        stop("if can mutate length is enabled max_length can't be NULL")
      }
      population[[i]] = mutateTree(population[[i]]$genes, percentage, TRUE, max_length)
    }else{
      population[[i]] = mutateTree(population[[i]]$genes, percentage)
    }
  }
  return(population)
}

init2 <- function(initial_solution, pop_size, mutation_percentage, max_genes){ # only used when an initial solution is provided
  # Create a population seeded with an initial solution, then mutate copies
  population = NULL
  population[[1]] = initial_solution
  write.table(2, file = "out.txt", sep = "\n",
              row.names = TRUE, col.names = NA)
  for(i in 2:pop_size){
    population[[i]] = mutateTree(population[[1]]$genes, mutation_percentage, TRUE, max_genes)
  }
  return(population)
}