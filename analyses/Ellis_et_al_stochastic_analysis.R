#### Ellis et al stochastic analysis

rm(list = ls() )
library(MASS)
setwd("~/Documents/Courses/Matrix Population Models/final_project/analyses/")
source("deterministic_analysis_tools.R")
source("stochastic_analysis_tools.R")

raw_dat = read.table("Ellis_et_al_Transition_Matrices.txt", header = TRUE, sep = "\t")

ordered_dat = raw_dat[ order(raw_dat$SPP, raw_dat$POP, raw_dat$YR), ]
names(ordered_dat)

population = data.frame( ordered_dat[, 1:3])

all_matrices = read_all_matrices( ordered_dat, 4)
all_matrices

Time = 5000

all_stochastic_results = iterate_stochastic_analysis(all_matrices, population, index_by = 2, Time)


results_df = compile_stochastic_results( all_stochastic_results, population)
results_df2 = merge(unique(population[, c(1,2)]), results_df,  by.all = POP)
results_df2
hist(results_df2$num_years)
hist(results_df2$sum_SE_mu)
hist(results_df2$sum_FE_mu)
hist(results_df2$sum_SE_var)
hist(results_df2$sum_FE_var)

write.csv(results_df2, "stochastic_analysis_results.csv", row.names = F)

