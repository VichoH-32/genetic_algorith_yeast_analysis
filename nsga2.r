## NSGA-II runner and plotting utilities for yeast gene selection ##
## Implements a domain-specific NSGA-II loop where individuals are trees
## built from sets of gene indices. Objectives used: RF-distance (score)
## and number-of-genes (length).

# path = "~/Documents/u2021/Tesis/final_runs/"
# path = "~/final_runs/"
if( !exists("path") ) path = "C:/Users/Vichi/Documents/U2021/"
# if( !exists("path") ) path = "~/Documents/U2021/"
library("phangorn")
library("MaOEA")
library("nsga2R")
library("emoa")
library("ggplot2")
library("BBmisc")
library("readxl")
library("MaOEA")
library("ecr")
print(1)
source(paste0(path,"code_base/base.r"))
print(2)
source(paste0(path,"code_base/mutation.r"))
source(paste0(path,"code_base/crossover.r"))

nsga2Yeast <- function(generations, pop_size, mutation_percentage, max_genes=100, initial_solution=NULL, plot=FALSE, images_dir, one.image=FALSE, fronts=NULL){
  # Main NSGA-II loop.
  # Inputs: generations, population size, mutation rate, max genes per solution
  # initial_solution: optional single tree to seed the population
  sols = NULL
  pareto_front_history = NULL
  write.table(1, file = "out.txt", sep = "\n",
              row.names = TRUE, col.names = NA)

  # Initialize parents: either random or seeded using init2
  parents = if(is.null(initial_solution)) init(pop_size, max_genes) else init2(initial_solution, pop_size, mutation_percentage, max_genes)
  # initialize offspring (will be replaced during evolution)
  offspring = init(pop_size, max_genes)

  for(gen_counter in 1:generations){
    #print(paste0("gen ", gen_counter))
    
    # Combine parents and offspring to form the population for selection
    population = c(parents, offspring)
    
    # prepare data to use in fastNonDominatedSorting
    datos = NULL
    for(i in 1:length(population)){
      datos = rbind(datos, c(population[[i]]$score, length(population[[i]]$genes), i, gen_counter))
    }
    colnames(datos) = c("score", "length", "position", "generation")

    # Rank solutions using nondominated sorting on the two objectives
    ranking = fastNonDominatedSorting(datos[,c(1,2)])

    # add ranking to each solution
    ranking_tmp = rep(0,nrow(datos))
    for (a in 1:length(ranking)){
      for (b in 1:length(ranking[[a]])){
        ranking_tmp[ranking[[a]][b]] = a
      }
    }

    datos = cbind(datos,ranking=ranking_tmp)
    # sort by ranking (front number) ascending
    datos = datos[order(datos[,"ranking"]),]
    
    #print(datos)
    
    # Calculate crowding distance per front so we can prefer diverse solutions
    crowd_dist = rep(Inf,nrow(datos))
    for (a in 1:length(unique(datos[,ncol(datos)]))){
      indices=which(datos[,ncol(datos)]==a)
      if (length(indices)>2){
        # emoa::crowding_distance expects points as rows; t(...) used to match signature
        crowd_dist[indices] = emoa::crowding_distance(t(datos[indices,c(1,2)]))
      }
    }
    datos = cbind(datos,crowd_dist=crowd_dist)

    # save pareto front info to history for plotting / HV computation later
    if(is.null(pareto_front_history)){
      #assign("pareto_front_history", datos, envir = .GlobalEnv)
      pareto_front_history = datos
    }else{
      #assign("pareto_front_history", rbind(pareto_front_history, datos), envir = .GlobalEnv)
      pareto_front_history = rbind(pareto_front_history, datos)
    }

    if(gen_counter == generations){
      # At the last generation build the list of best (rank 1) solutions
      df = data.frame(datos)
      best = df[df$ranking == 1,]
      for(i in 1:nrow(best)){
        sols[[i]] = population[[best[i, "position"]]]
	      sols[[i]]$ranking = best[i, "ranking"]
      }
    }else{
      # Order by ranking (asc) and crowding distance (desc) to pick parents
      datos = datos[order(datos[, "ranking"], -datos[,"crowd_dist"]),]

      # select the top-N positions to become the next generation parents
      n_pos = datos[1:pop_size,"position"]
      #print(n_pos)
      
      parents = NULL
      for(i in 1:pop_size){
        parents[[i]] = population[[n_pos[i]]]
      }

      # Create offspring by crossover + mutation
      if(length(offspring) > 1){
        offspring = crossover(parents, max_genes)
        offspring = mutation(offspring, mutation_percentage)
      }
    }
  }
  if(plot){
    generation.plot(pareto_front_history, images_dir, max_genes, one.image, fronts)
  }
  return(sols)
}

generation.plot.one <- function(data, images_dir, carpeta){
  # Plot the evolution of the Pareto front (single combined plot)
  df = data[data$ranking == 1,]
  print(df)
  g <- ggplot(df, aes(x=norm_score, y=norm_length, group = generation, color = generation)) + 
    geom_line(show.legend = FALSE) +
    geom_point(show.legend = FALSE) +
    xlab("Puntaje") + ylab("Largo") +
    expand_limits(x = 0, y = 0) + 
    scale_x_continuous(expand = c(0, 0), limits = c(0,1)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
    theme_bw()
  plot(g)
  ggsave(paste0(images_dir, "/SET", carpeta , "/nsga2_single", ".png"))
}

generation.plot.multi <- function(data, images_dir, carpeta){
  # Save one image per generation showing the Pareto front of rank 1
  for( i in 1:length(unique(data[, "generation"]))){
    df = data[data$generation == i & data$ranking == 1,]
    print(df)
    g <- ggplot(df, aes(x=norm_score, y=norm_length)) + 
      geom_line(show.legend = FALSE) +
      geom_point(show.legend = FALSE) +
      xlab("Puntaje") + ylab("Largo") +
      expand_limits(x = 0, y = 0) + 
      scale_x_continuous(expand = c(0, 0), limits = c(0,1)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
      theme_bw()
    plot(g)
    ggsave(paste0(images_dir, "/SET", carpeta , "/nsga2_", i, ".png"))
  }
}

generation.hv <- function(fronts, values, ref = c(1,1)){
  # offset <- 1
  # to <- length(values) 
  # for(i in (offset + 1):to) {
  #   prevFront = data.matrix(fronts[fronts$generation==values[i-offset], c("norm_score", "norm_length")])
  #   actualFront = data.matrix(fronts[fronts$generation == values[i], c("norm_score", "norm_length")])
  #   dist = HV( prevFront, actualFront)
  #   print(paste0("Distancia desde frontera ", values[i-offset], " a frontera ", values[i], ":"))
  #   print(dist)
  # }
  # Compute hypervolume for selected generations (values is a vector of generation numbers)
  for( i in values){
      print(paste0("HV frontera ", i, ":"))
      print(GetHypervolume( t(fronts[fronts$generation==i, c("norm_score", "norm_length")]), reference = ref ))
  }
  # for( i in 1:max(fronts[,"generation"])){
  #   print(paste0("HV frontera ", i, ":"))
  #   print(GetHypervolume( t(fronts[fronts$generation==values[i], c("norm_score", "norm_length")]), reference = ref ))
  # }
}

generation.plot <- function(current_data, images_dir, max_genes, one.image=FALSE, fronts=NULL){
  # Normalize objectives and create plot directory
  data = data.frame(current_data)
  
  objs = cbind(score=sapply(data[,"score"], as.numeric),length=sapply(data[,"length"], as.numeric))
  datos_norm = Normalize(t(objs))
  data = cbind(data, norm_score=datos_norm$normalizedObjective[1,], norm_length=datos_norm$normalizedObjective[2,])
  
  carpeta = length(list.files(images_dir)) + 1
  dir.create(paste0(images_dir, "/SET", carpeta))
  generation.plot.multi(data, images_dir, carpeta)
  if(one.image){generation.plot.one(data, images_dir, carpeta)}
  if(!is.null(fronts)){generation.hv(data, fronts)}
}

fronteras = nsga2Yeast(19, 14, 0.6, 22, NULL, TRUE, "C:/Users/Vichi/Documents/U2021/Images", one.image=TRUE, fronts=c(1,7,14,19))
# fronteras = nsga2Yeast(2, 5, 0.6, 10, NULL, TRUE, "C:/Users/Vichi/Documents/U2021/Images", one.image=TRUE, fronts=c(1,2))