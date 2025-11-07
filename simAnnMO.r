if( !exists("path") ) path = "~/U2021/"
# path = "~/final_runs/"
library("MaOEA")
library("ggplot2")
library("BBmisc")
print("a")
source(paste0(path,"code_base/base.r"))
print("b")
source(paste0(path,"code_base/mutation.r"))

checkDominance <- function(prev_sol, new_sol){ # only 2 objectives, positive: sol1 dominates, negative:sol2 dominates, 0:none dominates
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
  
  actual_solution = if(is.null(initial_solution)) buildTree(getRandomGenes(sample(1:max_genes, 1))) else initial_solution
  
  #generation_history = NULL
  
  for (i in 1:external_loops){
    for (j in 1:internal_loops){
      
      # mutate tree
      new_solution = mutateTree(actual_solution$genes, mutation_percentage, TRUE, max_genes)
      
      #Evaluación
      total_score = checkDominance(actual_solution, new_solution)
      
      if(total_score < 0){#Si la new_solution domina a actual_solution
        actual_solution = new_solution
      }else if(total_score > 0){#Si actual_solution domina a new_solution -> se calcula temperatura
        deltaE = new_solution$score - actual_solution$score
        if (checkTemperature(deltaE)){actual_solution = new_solution}
      }else{ # = 0, ninguna domina -> se calcula hipervolumen
        max_length = if(max_genes<length(actual_solution$genes)) length(actual_solution$genes) else max_genes
        sol1_len = length(actual_solution$genes) / max_length
        sol2_len = length(new_solution$genes) / max_length
        
        sols_data = cbind(c(actual_solution$score, new_solution$score), c(length(actual_solution$genes), length(new_solution$genes)), c(sol1_len, sol2_len))
        colnames(sols_data) = c("score", "length", "norm_length")
        hyp_point=GetHVContribution(t(sols_data[, c(1,3)]),reference = c(2,2))
        print(hyp_point)
        if( hyp_point[2] > hyp_point[1]){ # Si new_solution contribuye mas que actual_solution al HV
          actual_solution =  new_solution
        }else if(checkTemperature(hyp_point[2] - hyp_point[1])){ # se checkea temperatura
          actual_solution =  new_solution
        }
      }
      
      #plot part
      
    }
    T = T * alpha
  }
  return(actual_solution)
}

#init_sol = buildTreeFromXSLX("Resumen de Resultados.xlsx", 2)
init_sol = NULL

#external_loops, internal_loops, mutation_percentage, T, alpha, max_genes=100, initial_solution=NULL, plot=FALSE, images_dir
sol = simulatedAnnealingMO(external_loops=2, internal_loops=2, mutation_percentage=0.6, T=1000, alpha=0.8, initial_solution=init_sol)
