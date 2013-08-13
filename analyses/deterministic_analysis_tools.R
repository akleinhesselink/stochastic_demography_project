#### Deterministic matrix analysis functions ####
#### Adapted from code from Dave Koons' population biology class #### 

library(MASS)

calc_mean_matrix = function( list_of_matrices ){ 
  mat_dims = dim(list_of_matrices[[1]])
  
  total_matrix = matrix(0, mat_dims[1], mat_dims[2])
  for (i in list_of_matrices){
    total_matrix = i + total_matrix
  }
  mean_matrix = total_matrix/length(list_of_matrices)
  return (mean_matrix)
}

get_survival_sum = function( M ){ 
  elements =  M[-1, ] 
  survival_sum = sum(elements)
  return(survival_sum)
}

get_retrogression_sum = function( M ){ 
  elements = upper.tri(M, diag = F)  
  elements[1, ] = FALSE
  retro_sum = sum(M[elements])
  return(retro_sum)
}

get_subdiagonal_sum = function( M ){
  elements = lower.tri(M, diag = TRUE)
  sub_diag_sum = sum(M[elements])
  return(sub_diag_sum)  
}

get_fecundity_sum = function( M ){ 
  elements = M[1, ]
  fecundity_sum = sum(elements)
  return(fecundity_sum)
}

get_dominant_eigenvalues = function( A ){
  D<-dim(A)
  rows<-D[1]
  cols<-D[2]
  eig = eigen(A)                    # eigenvalues of A
  EigVecs = eig$vectors             # eigenvectors of A
  Lambdas = Re(eig$values)          # real number components of eigenvalues
  Lambda = max(Lambdas)             # dominant eigenvalue...long-term geometric rate of population growth
  pos = which.max(Lambdas)          # finding the position of the dominant eigenvalue amongst all of them
  w = Re(eig$vectors[1:rows,pos])   # its associated right eigenvector
  V = Conj(ginv(EigVecs))           # left eigenvector
  v = Re(t(t(V[pos,])))             # dominant left eigenvector

  out = list(Lambda, w, v)
  names(out) = c('Lambda', 'w', 'v')
  return(out)
}
 

get_sad = function ( A ){ 
  with( get_dominant_eigenvalues( A ), {
    sad = w/(sum(w))
    sad = round(sad, 3)
    return ( sad ) 
  })
}
  

get_rv = function ( A ){
  with( get_dominant_eigenvalues ( A ), { 
    rv = v/(sum(v))
    rv = round(rv,3)
    return( rv )
  })
}


get_life_expectancy = function( A ){
  Lambda = get_dominant_eigenvalues( A )[[1]]
  D<-dim(A)
  rows<-D[1]
  cols<-D[2]
  TRAN = rbind(0,A[2:rows,])
  FECU = rbind(A[1,],matrix(0,(rows-1),cols))  # this separates A into fertility (F) and transition matrices (T)
  n = c(rows)
  I = matrix(0,nrow=n,ncol=n)
  I[row(I)==col(I)] = 1             # Identity matrix
  N = ginv(I-TRAN)                  # something called the fundamental matrix
  R = FECU%*%N
  eigR = eigen(R)
  Rreal = Re(eigR$values)
  R0 = max(Rreal)                   # expected lifetime reproductive output (per generation growth rate)
  GT = log(R0)/log(Lambda)          # Time required to achieve R0, a measure of generation time
  D = ginv((I-TRAN)%*%(I-TRAN))
  mu1 = (t(D[,1])%*%FECU[1,])/R0       # Mean age of the parents for a cohort of newborns, another measure of generation time
  LE = colSums(N)                   # the mean life expectancy from the beginning of each stage
  
  out = list(R0, GT, mu1, LE)
  names(out) = c('R0', 'GT', 'mu1', 'LE')
  return (out) 
}


get_sensitivities = function ( A ) { 
  with( get_dominant_eigenvalues ( A ), { 
    senmat = v%*%t(w)                 # raw sensitivity matrix
    senmat[A==0] = 0                  # puts 0s in locations where vital rate does not exist
    return ( senmat )
  })
}


get_elasticities = function ( A ) { 
  with( get_dominant_eigenvalues ( A ), { 
    # Elasticity analysis of Lambda to proportional changes in vital rates
    emat = A/Lambda*(v%*%t(w))
    return ( emat )
  })
}

convert_mx <- function(mx_str){
  #### reads and converts matrix data from Ellis et al. 2012
  #### source code at: http://esapubs.org/archive/ecol/E093/083/metadata.htm
  mx_str<-substr(mx_str, 2, nchar(mx_str)-1)
  mx_str<-gsub('; ', ',', mx_str)
  mx_str<-gsub(' ', ',', mx_str)
  mx_str<-paste('c(',mx_str,')',sep='')
  mx <- eval(parse(text=mx_str))
  mx<-matrix(mx, nrow=sqrt(length(mx)), byrow=T)
  return(mx) 
}

read_all_matrices = function ( raw_data, trmx_col){
  all_matrices = list(NA)
  for (row in 1:nrow(raw_data)){
    all_matrices[[row]] = (convert_mx(as.character(raw_data[row, trmx_col])))
  }  
  return (all_matrices)
}

iterate_calc_mean_matrix = function( list_of_matrices, population_df, index_by ){ 
  mean_matrix_list = list(NA)
  row = 1
  for (p in levels(population_df[, index_by])){ 
    index = which(population_df[, index_by] == p )
    mean_matrix_list[[row]] = calc_mean_matrix(list_of_matrices[index]) 
    row = row + 1
  }
  return(mean_matrix_list)  
}


iterate_deterministic_analysis = function( list_of_matrices, population_df){ 
  row = 1
  deterministic = data.frame( unique(population_df[, 1:2]), Lambda = NA, R0 = NA, GT= NA, mu1 = NA, LE2 = NA)
  for ( mat in list_of_matrices ){ 
    deterministic[row, ]$Lambda = as.numeric(get_dominant_eigenvalues(mat)[1])
    deterministic[row, ]$R0 = get_life_expectancy(mat)$R0
    deterministic[row, ]$GT = get_life_expectancy(mat)$GT
    deterministic[row, ]$mu1 = as.numeric(get_life_expectancy(mat)$mu1)
    deterministic[row, ]$LE2 = as.numeric(get_life_expectancy(mat)$LE[2])
    row = row + 1
  }
  return(deterministic)  
}

iterate_deterministic_EMAT = function( list_of_matrices, pop_names){
  i = 1
  output = list(NA)
  for ( mat in list_of_matrices) { 
    output[[i]] = get_elasticities(mat)
    i = i + 1
  }
  names(output) = pop_names
  return(output)
}


