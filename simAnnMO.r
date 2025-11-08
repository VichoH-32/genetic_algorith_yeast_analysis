## Simulated annealing (multi-objective) helper ##
## This file contains a simplistic multi-objective simulated annealing implementation
## that operates on tree solutions (same solution representation as the GA code).
## The objective of this implementation is to have a comparation point with the NSGA-II code.

if( !exists("path") ) path = "~/U2021/"
# path = "~/final_runs/"
library("MaOEA")
library("ggplot2")
library("BBmisc")
#print("a")
source(paste0(path,"code_base/base.r"))
#print("b")
source(paste0(path,"code_base/mutation.r"))

checkDominance <- function(prev_sol, new_sol){ # only 2 objectives
  # returns: positive -> prev_sol dominates (old is better), negative -> new_sol dominates (new is better), 0 -> non-dominated
  total_score = 0
  total_score = if(prev_sol$score <= new_sol$score) (total_score + 1) else (total_score - 1)
  total_score = if(length(prev_sol$genes) <= length(new_sol$genes)) (total_score + 1) else (total_score - 1)
  return(total_score)
}

checkTemperature <- function(deltaE){
  prob = exp((-deltaE/T))
  prob_acep = runif(1)
  return(prob_acep < prob)
}
#simulatedAnnealingMO <- function(pop_size, external_loops, mutation_percentage, initial_solution=NULL, plot=FALSE, images_dir){
simulatedAnnealingMO <- function(external_loops, internal_loops, mutation_percentage, T, alpha, max_genes=100, initial_solution=NULL, plot=FALSE, images_dir){
  # external_loops: number of outer temperature reductions
  # internal_loops: number of iterations per temperature
  # T, alpha: initial temperature and cooling factor (T <- T * alpha each outer loop)

  actual_solution = if(is.null(initial_solution)) buildTree(getRandomGenes(sample(1:max_genes, 1))) else initial_solution
  
  #generation_history = NULL
  
  for (i in 1:external_loops){
    for (j in 1:internal_loops){
      # propose a neighbor by mutating the current solution
      new_solution = mutateTree(actual_solution$genes, mutation_percentage, TRUE, max_genes)

      # Determine dominance relationship
      total_score = checkDominance(actual_solution, new_solution)

      if(total_score < 0){
        # new_solution dominates -> accept
        actual_solution = new_solution
      }else if(total_score > 0){
        # actual_solution dominates -> accept with temperature probability check
        deltaE = new_solution$score - actual_solution$score
        if (checkTemperature(deltaE)){actual_solution = new_solution}
      }else{ # neither dominates -> use hypervolume contribution as tie-breaker
        max_length = if(max_genes<length(actual_solution$genes)) length(actual_solution$genes) else max_genes
        sol1_len = length(actual_solution$genes) / max_length
        sol2_len = length(new_solution$genes) / max_length

        sols_data = cbind(c(actual_solution$score, new_solution$score), c(length(actual_solution$genes), length(new_solution$genes)), c(sol1_len, sol2_len))
        colnames(sols_data) = c("score", "length", "norm_length")
        hyp_point=GetHVContribution(t(sols_data[, c(1,3)]),reference = c(2,2))
        print(hyp_point)
        if( hyp_point[2] > hyp_point[1]){ # new solution contributes more to HV -> accept
          actual_solution =  new_solution
        }else if(checkTemperature(hyp_point[2] - hyp_point[1])){ # temperature probability check
          actual_solution =  new_solution
        }
      }
    }
    # cool down temperature
    T = T * alpha
  }
  return(actual_solution)
}

#init_sol = buildTreeFromXSLX("Resumen de Resultados.xlsx", 2)
init_sol = NULL

# Example usage (commented):
#external_loops, internal_loops, mutation_percentage, T, alpha, max_genes=100, initial_solution=NULL, plot=FALSE, images_dir
sol = simulatedAnnealingMO(external_loops=2, internal_loops=2, mutation_percentage=0.6, T=1000, alpha=0.8, initial_solution=init_sol)
