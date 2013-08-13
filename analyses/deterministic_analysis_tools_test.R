##################################################################
# Adapted from Dave Koons Class on Matrix Models 
##################################################################

rm(list=ls(all=TRUE))

setwd("~/Documents/Courses/Matrix Population Models/final_project")

source('deterministic_analysis_tools.R')

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

list_of_matrices = list(A_1, A_2, A_3, A_4)

#### The expected values are shown in comments below: 

A_test = calc_mean_matrix ( list_of_matrices )
A_test
#### 0.0, 0.5250
#### 0.6, 0.6875

get_dominant_eigenvalues(A_test)
#### Lambda = 1.001902  
#### w = -0.4641417, -0.8857610 
#### v =       
#### -0.5146103
#### -0.8593152

get_rv(A_test)
#### 0.375 
#### 0.625

get_sad(A_test)
#### 0.344 0.656

get_life_expectancy(A_test)
#### R0 = 1.008
#### GT = 4.193334
#### LE = 2.92, 3.20
#### mu1 = 4.2

get_sensitivities(A_test)
####     [,1]  [,2]
####[1,] 0.000 0.456
####[2,] 0.399 0.761

get_elasticities(A_test)
####      [,1]  [,2]
####[1,] 0.000 0.239
####[2,] 0.239 0.522



