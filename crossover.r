if( !exists("path") ) path = "~/Documents/U2021/"
source(paste0(path,"code_base/base.r"))

### Global variables ###

### Functions ###

crossSolutions <- function(sol1, sol2, max_genes){
  if(length(sol1) == 1){sol1 = append( sol1, getRandomGenes(1))}
  if(length(sol2) == 1){sol2 = append( sol2, getRandomGenes(1))}
  
  half1 = sample(1:(length(sol1)-1), 1)
  half2 = sample(1:(length(sol2)-1), 1)
  
  list1 = unique(append( sol1[1:half1], sol2[(half2+1):length(sol2)]))
  list2 = unique(append( sol2[1:half2], sol1[(half1+1):length(sol1)]))
  
  if(length(list1) > max_genes ){list1 = unique(sample(list1,max_genes))}
  if(length(list2) > max_genes ){list2 = unique(sample(list2,max_genes))}
  
  newSol1 = buildTree(list1)
  newSol2 = buildTree(list2)
  
  if(newSol1$score < newSol2$score){return(newSol1)}else{return(newSol2)}
}

crossover <- function(population, max_genes=0){
  newPopulation = NULL
  for(i in 1:length(population)){
    randomSols = sample(1:length(population), 2)
    newPopulation[[i]] = if(max_genes > 0) crossSolutions(population[[randomSols[1]]]$genes, population[[randomSols[2]]]$genes, max_genes) else crossSolutions(population[[randomSols[1]]]$genes, population[[randomSols[2]]]$genes)
    #newPopulation[[i]] = crossSolutions(population[[randomSols[1]]]$genes, population[[randomSols[2]]]$genes)
  }
  return(newPopulation)
}