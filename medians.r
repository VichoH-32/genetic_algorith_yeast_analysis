## Utilities to aggregate and compare results stored as .RData files ##
## Loads saved `res` objects produced by the runner, builds a table of scores
## and computes quick statistics (medians, standard deviations) per folder.
path = "C:\\Users\\Vichi\\Documents\\U2021\\final_runs_results\\"
folders = c("nsga2_short", "testNew")
df = NULL

for( folder in folders){
  files = list.files(paste0(path,folder))
  for(file in files){
    if(grepl("\\.RData", file)){
      load(paste0(path, folder, "\\", file))
      for(i in 1:length(res)){
        if(is.null(df)){
          # df = data.frame(score=res[[i]]$score, length=length(res[[i]]$genes), in_file_pos=i, file_name=file, algorithm=algorithm)
          df = data.frame(score=res[[i]]$score, length=length(res[[i]]$genes), in_file_pos=i, folder=folder, file=file)
        }else{
          # df = rbind(df, c(res[[i]]$score, length(res[[i]]$genes), i, file, algorithm))
          df = rbind(df, c(res[[i]]$score, length(res[[i]]$genes), i, folder, file))
        }
      }
    }
  }
}

print(df)

prev = sapply(df[ df$folder == "nsga2_short","score"], as.numeric)
new = sapply(df[ df$folder == "testNew","score"], as.numeric)

median(prev)
median(new)

sd(prev)
sd(new)