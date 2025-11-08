## Crossover operators for genetic algorithm ##
## This file provides two functions:
## - crossSolutions: produce two children (and pick the better) from two parents
## - crossover: apply crossSolutions repeatedly to build a new offspring population

if( !exists("path") ) path = "~/Documents/U2021/"
source(paste0(path,"code_base/base.r"))

### Functions ###

crossSolutions <- function(sol1, sol2, max_genes){
  # sol1 / sol2 are vectors of gene indices. If a parent has only one gene
  # we append one random gene to avoid degenerate crossover.
  if(length(sol1) == 1){sol1 = append( sol1, getRandomGenes(1))}
  if(length(sol2) == 1){sol2 = append( sol2, getRandomGenes(1))}
  
  # pick a random cut-point for each parent (simple one-point crossover)
  half1 = sample(1:(length(sol1)-1), 1)
  half2 = sample(1:(length(sol2)-1), 1)
  
  # build children by concatenating parts from each parent, then unique-ify
  list1 = unique(append( sol1[1:half1], sol2[(half2+1):length(sol2)]))
  list2 = unique(append( sol2[1:half2], sol1[(half1+1):length(sol1)]))
  
  # if child exceeds allowed max_genes, randomly sample down to max_genes
  if(length(list1) > max_genes ){list1 = unique(sample(list1,max_genes))}
  if(length(list2) > max_genes ){list2 = unique(sample(list2,max_genes))}
  
  # build trees for both children and return the better (lower score)
  newSol1 = buildTree(list1)
  newSol2 = buildTree(list2)
  if(newSol1$score < newSol2$score){
    return(newSol1)
  }else{
    return(newSol2)
  }
}

crossover <- function(population, max_genes=0){
  # Given a population (list of tree objects), produce a new population by
  # crossing randomly selected parents.
  newPopulation = NULL
  for(i in 1:length(population)){
    randomSols = sample(1:length(population), 2)
    if(max_genes > 0) {
      newPopulation[[i]] = crossSolutions(population[[randomSols[1]]]$genes, population[[randomSols[2]]]$genes, max_genes)
    } else {
      newPopulation[[i]] = crossSolutions(population[[randomSols[1]]]$genes, population[[randomSols[2]]]$genes)
    }
    #newPopulation[[i]] = crossSolutions(population[[randomSols[1]]]$genes, population[[randomSols[2]]]$genes)
  }
  return(newPopulation)
}