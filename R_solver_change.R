################ main solvers ##################
# classic lasso solver  
ggd_solver <- function(Y, mask, lambda, theta = NULL, num_it, thresh){
  if(is.null(theta)){
    theta <- Y * 0
  }
  norm_diff <- rep(0,num_it)
  
  for(it in 1:num_it){
    U <- Y * mask + theta * (1-mask)
    theta_new <- singular_soft_thresh(U,lambda)
    norm_diff[it] <- mean((theta_new - theta)^2)
    theta <- theta_new
    if(norm_diff[it] < thresh) break;
  }
  return(list(theta=theta, norm_diff = norm_diff))
}

# squared root lasso solver 
sqrt_ggd_solver <- function(Y, mask, lambda, theta = NULL, num_it, thresh, sigma = 1 ){
  if(is.null(theta)){
    theta <- Y * 0
  }
  norm_diff <- rep(0,num_it)
  
  for(it in 1:num_it){
    U <- Y * mask + theta * (1-mask)
    theta_new <- singular_soft_thresh(U,lambda * sigma)
    norm_diff[it] <- mean((theta_new - theta)^2)
    sigma_new = norm_obs(Y - theta_new * mask)
    sigma = sigma_new
    theta = theta_new
    if(norm_diff[it] < thresh) break;
  }
  return(list(theta=theta, norm_diff = norm_diff,sigma = sigma))
}
################ end of main solvers ##################

############## little functions ######################
# singular soft threshold for the matrice 
singular_soft_thresh <- function(U, lambda){
  junk <- svd(U)
  new_singular_vals <- pmax(junk$d - lambda, 0)
  new_U <- junk$u %*% diag(new_singular_vals) %*% t(junk$v)
  return(new_U)
}

# Mask function 
## this is for later use, Here the mask is generated with the data, so no need.   
## mask function - indicator 0/1 matrix of observed/unobseved entries 
id_obs = function(x){
  vec_x = as.vector(x)
  vec_x[which(vec_x != 0)] <- 1
  
  return(matrix(vec_x,nrow = nrow(x), ncol=ncol(x)))
}

# generating random rank-specified sparse matrix
gen_data <- function(n, rank, percent_missing, sigma){
  mu_init <- matrix(rnorm(n*n,mean = 3), ncol = n)
  junk <- svd(mu_init)
  if(rank != 1) mu = junk$u[,1:rank] %*% diag(junk$d[1:rank]) %*% t(junk$v[,1:rank])
  if(rank == 1 ) mu = junk$u[,1] %*% t(junk$v[,1]) * junk$d[1] 
  Y_init <- mu + rnorm(n*n,sd = sigma)
  mask <- matrix(rbinom(n*n, 1, 1-percent_missing), ncol = n)
  Y_obs <- Y_init * mask
  
  return(list(mask = mask, Y = Y_obs, mu = mu, SNR = norm(mu, type = "F")/(n*sigma)))
  #that is the same as calculating SNR <- norm(gen$mu, type = "F")/ norm(noise, type = "F")
}


# Frobenius norm for the non-zeros (observed) entries 
norm_obs<-function(x){
  n_obs <- length(which(x!=0))
  n <- nrow(x) * ncol(x)
  return(sqrt(norm(x,type = "F")^2 * n / n_obs ))
}

#In my implementation I wrote two little helper functions (to determine the foldids and the lambda_sequence
#if none are specified)
create.lambda <- function(data){
  return(seq(0,max(svd(data)$d),by = 1))
}
create.folds <- function(data, mask, nfold){
  length <- nrow(data)^2
  fold_vector <- sample(1:nfold, length, replace = TRUE)
  fold_matrix <- matrix(fold_vector, nrow = nrow(data)) * mask
  return(fold_matrix)
}
############## end of little functions ######################

########## lasso ########
lasso_mc <- function(data,mask,true,theta = NULL, lambda_sequence = NULL, thresh = 0.01){
  if(is.null(lambda_sequence)) lambda_sequence <- create.lambda(data)
  error.by.lambda <- numeric(length(lambda_sequence))
  unobserv.error.by.lambda <-numeric()
  observ.error.by.lambda <- numeric()
  # warm start for the first lambda 
  tst=theta
  for(nlambda in 1:length(lambda_sequence)){
    results <- ggd_solver(Y = data,theta = tst, mask = mask, lambda = lambda_sequence[nlambda],num_it = 1000,thresh = thresh )
    tst <- results$theta  
    unobserv.error.by.lambda[nlambda] <- norm(tst * (1-mask) - true * (1-mask),type = "F")
    observ.error.by.lambda[nlambda] <- norm(tst * mask - true * mask,type = "F")
    error.by.lambda[nlambda] <- norm(tst - true,type = "F") 
  }
  return(list(lambdas = lambda_sequence, 
              unobserv.errors = unobserv.error.by.lambda,
              observ.errors = observ.error.by.lambda,
              errors = error.by.lambda))
} 

########## end of lasso #######

######### square root lasso ######

sqrlasso_mc <- function(data,mask,true,theta = NULL, lambda_sequence = NULL, thresh = 0.01){
  if(is.null(lambda_sequence)) lambda_sequence <- create.lambda(data)
  error.by.lambda <- numeric()
  unobserv.error.by.lambda <-numeric()
  observ.error.by.lambda <- numeric()
  # warm start for the first lambda 
  tst=theta
  for(nlambda in 1:length(lambda_sequence)){
    results <- sqrt_ggd_solver(Y = data,theta = tst, mask = mask, lambda = lambda_sequence[nlambda],num_it = 1000,thresh = thresh )
    tst <- results$theta  
    unobserv.error.by.lambda[nlambda] <- norm(tst * (1-mask) - true * (1-mask),type = "F")
    observ.error.by.lambda[nlambda] <- norm(tst * mask - true * mask,type = "F")
    error.by.lambda[nlambda] <- norm(tst - true,type = "F") 
  }
  return(list(lambdas = lambda_sequence, 
              unobserv.errors = unobserv.error.by.lambda,
              observ.errors = observ.error.by.lambda,
              errors = error.by.lambda))
} 
######### end of square root lasso #######

############ cross validation ##############################
# cross valadation for classic lasso 
cv_matrixcmplt_ggd <- function(data,mask,true,theta = NULL, foldid = NULL, lambda_sequence = NULL, nfold = 10, thresh = 0.01){
  if(is.null(lambda_sequence)) lambda_sequence <- create.lambda(data)
  if(is.null(foldid))          foldid <- create.folds(data, mask, nfold)
  
  error.by.fold <- numeric(nfold)
  error.by.lambda <- numeric(length(lambda_sequence))
  
  # warm start for the fold 1
  for(nlambda in 1:length(lambda_sequence)){
    mask_tst = foldid
    mask_tst[ (mask_tst == 1) ] = 0
    mask_tst[ (mask_tst != 0) ] = 1
    results <- ggd_solver(Y = data * mask_tst,theta = theta, mask = mask_tst,lambda = lambda_sequence[nlambda],num_it = 1000,thresh = thresh )
    tst <- results$theta  
    error.by.fold[1] <- norm(tst * (1-mask_tst) - true * (1-mask_tst), type = "F") 
    ## using the initial beta as the warm start (for the lefting folds other than fold 1 )
    for(fold in 2:nfold){
      mask_tst = foldid
      mask_tst[ (mask_tst == fold) ] = 0
      mask_tst[ (mask_tst != 0) ] = 1
      results <- ggd_solver(Y = data * mask_tst,mask = mask_tst, theta = tst,lambda = lambda_sequence[nlambda],num_it = 1000,thresh = thresh )
      tst <- results$theta 
      error.by.fold[fold] <- norm(tst * (1-mask_tst) - true * (1-mask_tst),type = "F")
    }
    error.by.lambda[nlambda] <- mean(error.by.fold) 
  }
  return(list(lambdas = lambda_sequence, errors = error.by.lambda ))
} 


# cross validation for square-root lasso 
cv_matrixcmplt_ggd_sqr <- function(data, mask,true, theta = NULL, foldid = NULL, lambda_sequence = NULL, nfold = 10, thresh = 0.01){
  if(is.null(lambda_sequence)) lambda_sequence <- create.lambda(data)
  if(is.null(foldid))          foldid <- create.folds(data, mask, nfold)
  
  error.by.fold <- numeric(nfold)
  error.by.lambda <- numeric(length(lambda_sequence))
  
  # warm start for the fold 1
  for(nlambda in 1:length(lambda_sequence)){
    mask_tst = foldid
    mask_tst[ (mask_tst == 1) ] = 0
    mask_tst[ (mask_tst != 0) ] = 1
    results <- sqrt_ggd_solver(Y = data * mask_tst,theta = theta, mask = mask_tst,lambda = lambda_sequence[nlambda],num_it = 1000,thresh = thresh )
    tst <- results$theta  
    error.by.fold[1] <- norm(tst * (1-mask_tst) - true * (1-mask_tst),type = "F")
    
    ## using the initial beta as the warm start (for the lefting folds other than fold 1 )
    for(fold in 2:nfold){
      mask_tst = foldid
      mask_tst[ (mask_tst == fold) ] = 0
      mask_tst[ (mask_tst != 0) ] = 1
      results <- sqrt_ggd_solver(Y = data * mask_tst, mask = mask_tst, theta = tst,lambda = lambda_sequence[nlambda],num_it = 1000,thresh = thresh )
      tst <- results$theta 
      error.by.fold[fold] <- norm(tst * (1-mask_tst) - true * (1-mask_tst),type = "F")
    }
    error.by.lambda[nlambda] <- mean(error.by.fold) 
  }
  return(list(lambdas = lambda_sequence, errors = error.by.lambda ))
} 

############ end of cross validation ###########################

