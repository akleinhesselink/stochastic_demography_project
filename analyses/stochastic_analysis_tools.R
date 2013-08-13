##################################################################
# Adapted from Dave Koons Class on Matrix Models 
##################################################################

##################################################################
# Stochastic analysis with observed transition matrices 

source('deterministic_analysis_tools.R')
find_env_vector = function( list_of_matrices, Time ){
  num_mat = length(list_of_matrices)
  env = sample(num_mat, Time, replace = T)
  return (env)  
}


find_rv_vector = function( list_of_matrices, Time, env ){ 
  matdims = dim(list_of_matrices[[1]])
  rows = matdims[1]
  cols = matdims[2]
  v = matrix(0,rows,(Time+1))
  vvec = (matrix(1,1,cols))/cols
  v[,(Time+1)]=t(vvec)
  for (time in Time:1) {
    A = list_of_matrices[[env[time]]]
    vvec=vvec%*%A
    vvec=vvec/sum(vvec)
    v[,time]=t(vvec)
  }
  out = list(vvec, v)
  names(out) = c('vvec', 'v')
  return(out)
}


find_lambda = function(list_of_matrices, Time, env ){ 
  matdims = dim(list_of_matrices[[1]])
  rows = matdims[1]
  cols = matdims[2]
  wvec = (matrix(1,rows,1))/rows
  w = matrix(0,rows,(Time+1))
  w[,1] = wvec
  lam = matrix(0,Time,1)
  
  for (time in 1:Time) {
    A = list_of_matrices[[env[time]]]
    wvec=A%*%wvec
    lam[time]=sum(wvec)           # per time step population growth rate
    wvec=wvec/lam[time]           # per time step measures of age distribution
    w[,(time+1)]=wvec
  }
   
  ln_lambda_s = mean(log(lam))      # the stochastic population growth rate
  Lambda_s = exp(ln_lambda_s)       # the stochastic population growth rate on the
  out = list(Lambda_s, lam, w)
  names(out) = c('Lambda_s', 'lam', 'w')
  return(out)
}

get_stochastic_eigenvalues = function( list_of_matrices, Time, env){
  rv_out =  find_rv_vector (list_of_matrices, Time, env)
  vvec = rv_out$vvec
  v = rv_out$v
  lambda_out =  find_lambda( list_of_matrices, Time, env)
  Lambda_s = lambda_out$Lambda_s
  lam = lambda_out$lam
  w = lambda_out$w  
  out = list(env, rv_out, vvec, v, Lambda_s, lam, w)
  names(out) = c('env', 'rv_out', 'vvec', 'v', 'Lambda_s', 'lam', 'w')
  return (out)
}


get_one_EMAT = function( A, lambda, v_next, w, w_next ){ 
  #### BASED on formula:  (v[,time+1] %*% t(w[,time]) * A) / as.numeric(lam[time] * t(v[,time+1]) %*% w[,time+1])          
  EMAT = (v_next%*%t(w)*A)/as.numeric(lambda*t(v_next)%*%w_next)
  return ( EMAT )  
}


get_stochastic_elasticity = function ( list_of_matrices, Time, env){ 
  
  A_mean = calc_mean_matrix(list_of_matrices)
  matdims = dim(A_mean)
  rows = matdims[1]
  cols = matdims[2]
  
  #Matrices to store elasticities for matrix elements
  EMAT = matrix(0,rows,cols)        #Stochastic elasticity
  EMAT_mu = matrix(0,rows,cols)     #Elasticity to mean
  EMAT_sig = matrix(0,rows, cols)    #Elasticity to variance
  
  #### Get eigenvalues for each time step
  stoch_eig =  get_stochastic_eigenvalues(list_of_matrices, Time, env)
  
  with(stoch_eig, {   
    #Loop through years, starting at year 1000 to discard transience, and calculate elasticities
    for (time in 1000:Time){
      A = list_of_matrices[[env[time]]]
      sigma = A - A_mean                    #The difference between this year and the average year
      v_next = v[ , (time + 1)]
      w_now = w[, time ]
      w_next = w[, (time + 1)]
      lambda = lam[ time ]      
      EMAT_mu = EMAT_mu + get_one_EMAT(A_mean, lambda, v_next, w_now, w_next)
      EMAT_sig = EMAT_sig + get_one_EMAT(sigma, lambda, v_next, w_now, w_next) 
      EMAT = EMAT + get_one_EMAT(A, lambda, v_next, w_now, w_next)
      }
  
  EMAT = EMAT/length(1000:Time)
  EMAT_mu = EMAT_mu/length(1000:Time)
  EMAT_sig = EMAT_sig/length(1000:Time)
  
  out = list(EMAT, EMAT_mu, EMAT_sig)
  names(out) = c('EMAT', 'EMAT_mu', 'EMAT_sig')
  return ( out )
  })
}
  
  
stochastic_matrix_analysis = function ( list_of_matrices, sequential = F, Time) {
  if(sequential == F){
    env = find_env_vector(list_of_matrices, Time)
  }
  else if(sequential == T) {
    env = rep(1:length(list_of_matrices), length.out = Time)
  }
  rv_out = find_rv_vector(list_of_matrices, Time, env)
  vvec = rv_out[[1]]
  v = rv_out[[2]]
  RV = rowMeans(v) 
  Lambda_out = find_lambda(list_of_matrices, Time, env)
  Lambda_s = Lambda_out[[1]]
  lam = Lambda_out[[2]]
  w = Lambda_out[[3]]
  SAD = rowMeans(w)
  EMAT_s = get_stochastic_elasticity(list_of_matrices, Time, env)
  output = list(Lambda_s, SAD, RV, EMAT_s$EMAT, EMAT_s$EMAT_mu, EMAT_s$EMAT_sig)
  names(output) = c('Lambda_s', 'SAD_s', 'RV_s', 'EMAT', 'EMAT_mu', 'EMAT_sig')
  return(output)
}


iterate_stochastic_analysis = function( list_of_matrices, population_df, Time, index_by){ 
  results_list = list(NA)
  i = 1
  for (p in levels(population_df[, index_by])){ 
    index = which(population_df[, index_by] == p)
    results_list[[i]] = stochastic_matrix_analysis(list_of_matrices[index], sequential = F, Time) 
    i = i + 1
  }
  return(results_list)
}


old_stochastic_elasticity = function( list_of_matrices, Time, env){ 
  matdims = dim(list_of_matrices[[1]])
  rows = matdims[1]
  cols = matdims[2]
  
  stoch_eig =  get_stochastic_eigenvalues(list_of_matrices, Time, env)
  with( stoch_eig, { 
    EMAT = matrix(0,rows,cols)  
    for (time in 1:Time) {
      A = list_of_matrices[[env[time]]]
      v_next = v[ , (time + 1)]
      w_now = w[, time ]
      w_next = w[, (time + 1)]
      lambda = lam[ time ]
      EMAT =  EMAT + get_one_EMAT( A, lambda, v_next, w_now, w_next)    
    }
    EMAT = EMAT/Time
    return (EMAT)
  })
}

compile_stochastic_results = function( all_stochastic_results, population){
  results_df = data.frame(NA)
  j = 1
  for (i in all_stochastic_results){ 
    results_df[j, 1] = POP_name = levels(population[, 2])[j]
    results_df[j, 2] = i$Lambda_s   
    results_df[j, 3] = sum(i$EMAT_mu)
    results_df[j, 4] = sum(i$EMAT_sig)
    results_df[j, 5] = get_fecundity_sum(i$EMAT_mu)
    results_df[j, 6] = get_fecundity_sum(i$EMAT_sig)
    results_df[j, 7] = get_survival_sum(i$EMAT_mu)
    results_df[j, 8] = get_survival_sum(i$EMAT_sig)
    results_df[j, 9] = length( which(population[, 2] == POP_name))
    j = j + 1
  }
  names(results_df) = c('POP', 'Lambda_s', 'sum_E_mu', 'sum_E_var', 'sum_FE_mu', 'sum_FE_var', 
                        'sum_SE_mu', 'sum_SE_var', 'num_years')
  return( results_df)
}
