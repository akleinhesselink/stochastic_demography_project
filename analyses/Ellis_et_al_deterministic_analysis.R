#### Ellis et al. Deterministic analysis 
rm(list = ls() )
library(MASS)
setwd("~/Documents/Courses/Matrix Population Models/final_project/analyses/")
source("deterministic_analysis_tools.R")

raw_dat = read.table("Ellis_et_al_Transition_Matrices.txt", header = TRUE, sep = "\t")

ordered_dat = raw_dat[ order(raw_dat$SPP, raw_dat$POP, raw_dat$YR), ]
names(ordered_dat)

population = data.frame( ordered_dat[, 1:3])
plot(population$YR, type = "l") #### check that years are in order

all_matrices = read_all_matrices( ordered_dat, 4)
all_matrices

#### Deterministic analysis for each population's mean matrix

#### Test that orders match up 
for (i in 1:nrow(population)){
  test = population$POP[i]
  print (as.character(test))
  index = which(population$POP == test)
  test_mats = read_all_matrices(ordered_dat[ index, ], 4)
  print (try ( test_mats[[1]] == all_matrices[[index[1]]]))
  print (try ( test == ordered_dat[ index , 2]))
}

mean_matrices = iterate_calc_mean_matrix(all_matrices, population, index_by = 2)
deterministic_df = iterate_deterministic_analysis(mean_matrices, population)
deterministic_df$LE2[deterministic_df$LE2 < 0] = NA

EMAT_d = iterate_deterministic_EMAT(mean_matrices, levels(population$POP))
j = 1
for ( m in EMAT_d) { 
  deterministic_df$sum_SE[j] = get_survival_sum(m)
  deterministic_df$sum_FE[j] = get_fecundity_sum(m)  
  j = j+1 
}

deterministic = as.data.frame(deterministic_df)
write.csv(deterministic, "deterministic_analysis_results.csv", row.names = F)

