if( !exists("path") ) path = "~/Documents/U2021/"
#path = "~/final_runs/"
library(foreach)
library(doParallel)

source(paste0(path,"code_base/used_packages.r"))
source(paste0(path,"code_base/base.r"))
source(paste0(path,"code_base/nsga2.r"))

runner <- function(i, initial_solution=NULL){
  write.table(0, file = "out.txt", sep = "\n",
              row.names = TRUE, col.names = NA)
  res = nsga2Yeast(generations=19,
                   pop_size=14,
                   mutation_percentage=0.6,
                   max_genes=22,
                   initial_solution=initial_solution
  )
  filename = paste0(path,"medoides/nsga2_medoides",i,".RData")
  save(res, file=filename)
  print(paste0("saved ", filename))
  return(res)
}

base_tree = buildTreeFromXSLX("Grupos y medoides.xlsx", 4)
# base_tree = NULL
if(!is.null(base_tree)){
  base_tree$og = TRUE
  print("Base Tree Built")
}

registerDoParallel(5)
foreach(i=1:31, .packages = used_packages) %dopar% {runner(i, base_tree)}
