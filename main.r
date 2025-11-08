## Main runner for NSGA-II algorithm for yeast analysis ##
## This file contains the script entry point used to launch many
## independent NSGA-II runs in parallel and save their results.

# Default working path used by the scripts; projects typically set this
# externally before sourcing this file. It points to where data/code live.
if( !exists("path") ) path = "~/Documents/U2021/"
#path = "~/final_runs/" # alternative example

library(foreach)    # parallel foreach construct
library(doParallel) # backend for foreach

# Source shared helper scripts (package list, base utilities, NSGA-II impl)
source(paste0(path,"code_base/used_packages.r"))
source(paste0(path,"code_base/base.r")) 
source(paste0(path,"code_base/nsga2.r"))

# runner: wrapper executed in parallel. i is just an identifier to name the
# output file. initial_solution can be passed to seed the population.
runner <- function(i, initial_solution=NULL){
  # quick progress stamp to a file (used for debugging/monitoring)
  write.table(0, file = "out.txt", sep = "\n",
              row.names = TRUE, col.names = NA)

  # call the algorithm with chosen hyperparameters (example values)
  res = nsga2Yeast(generations=19,
                   pop_size=14,
                   mutation_percentage=0.6,
                   max_genes=22,
                   initial_solution=initial_solution
  )

  # save result object to disk (filename depends on 'i')
  filename = paste0(path,"medoides/nsga2_medoides",i,".RData")
  save(res, file=filename)
  print(paste0("saved ", filename))
  return(res)
}

# Optionally build a base tree from an excel file that contains predefined sets
# buildTreeFromXSLX returns a tree object compatible with the rest of the code
base_tree = buildTreeFromXSLX("Grupos y medoides.xlsx", 4)
# base_tree = NULL
if(!is.null(base_tree)){
  base_tree$og = TRUE # mark it so callers can detect it came from file
  print("Base Tree Built")
}

# Register parallel backend: 5 worker processes in this example.
# Tune this number to the available CPU cores / memory.
registerDoParallel(5)

# Launch 31 parallel runs. .packages ensures each worker loads required packages
foreach(i=1:31, .packages = used_packages) %dopar% {runner(i, base_tree)}
