## List of packages required by the analysis scripts ##
## This vector is used by the parallel runner to ensure workers load
## the same package set. Note: some packages appear duplicated below.

used_packages = c("phangorn", "MaOEA", "nsga2R", "emoa", "ggplot2", "BBmisc", "readxl", "MaOEA", "readxl")
# If you need to install missing packages uncomment the line below.
# install.packages(unique(used_packages))