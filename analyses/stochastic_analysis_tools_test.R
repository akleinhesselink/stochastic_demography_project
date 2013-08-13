##################################################################
# Adapted from Dave Koons Class on Matrix Models 
##################################################################

rm(list=ls(all=TRUE))

setwd("~/Documents/Courses/Matrix Population Models/final_project")
source('deterministic_analysis_tools.R')
source('stochastic_analysis_tools.R')

# matrix for 1st year of study
A_1 = matrix(0,2,2)
A_1[1,2] = 1.00
A_1[2,1] = 0.65
A_1[2,2] = 0.75

# matrix for 2nd year of study
A_2 = matrix(0,2,2)
A_2[1,2] = 0.75
A_2[2,1] = 0.55
A_2[2,2] = 0.65

# matrix for 3rd year of study
A_3 = matrix(0,2,2)
A_3[1,2] = 0.25
A_3[2,1] = 0.45
A_3[2,2] = 0.50

# matrix for 4th year of study
A_4 = matrix(0,2,2)
A_4[1,2] = 0.10
A_4[2,1] = 0.75
A_4[2,2] = 0.85
A_m=(A_1+A_2+A_3+A_4)/4     # mean matrix

list_of_matrices = list(A_1, A_2, A_3, A_4)

Time = 100000  

env = find_env_vector (list_of_matrices, Time)
stoch_met = get_stochastic_eigenvalues( list_of_matrices, Time, env)

stoch_met$Lambda_s    #### Lambda_s
rowMeans(stoch_met$v) #### RV
rowMeans(stoch_met$w) #### SAD

EMAT = old_stochastic_elasticity(list_of_matrices, Time, env)
EMAT

#### Output from Dave Koon's original example as follows: 
#### Lambda_s = 0.97
#### SAD = 0.31, 0.69 
#### RV = 0.38, 0.62
#### EMAT = 0.00, 0.23 
####        0.23, 0.54
#### compare to output above

EMAT_s = get_stochastic_elasticity(list_of_matrices, Time, env)
EMAT_s

EMAT_s$EMAT - (EMAT_s$EMAT_mu + EMAT_s$EMAT_sig) #### should be close to zero 

all_out = stochastic_matrix_analysis(list_of_matrices, sequential = F, Time)
all_out


#### Compare with deterministic analsyis below
mean_matrix = calc_mean_matrix(list_of_matrices)
get_dominant_eigenvalues(mean_matrix) 
get_sad(mean_matrix)
get_rv(mean_matrix)
get_elasticities(mean_matrix)


