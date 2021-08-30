# cholesky2 <- function(matrix, toler=1e-5){
#
#  eps = 0
#  nonneg = 1
#
#  n <- ncol(matrix)
#
#  for (i in seq(n)){
#   if(matrix[i,i] > eps) eps <- matrix[i,i]
#   # copy upper right values to bottom left
#   for(j in seq(i+1, n)){
#    if(j <= n && j > i)
#     matrix[j,i] <- matrix[i,j]
#   }
#  }
#  if (eps == 0)
#   eps = toler # no positive diagonals!
#  else
#   eps = eps * toler
#
#  rank = 0
#
#  for (i in seq(n)) {
#   pivot = matrix[i, i]
#   if (!is.finite(pivot) || pivot < eps) {
#    matrix[i, i] = 0
#    if (pivot < -8*eps) nonneg= -1
#   } else  {
#    rank <- rank + 1
#    for(j in seq(i+1, n)){
#     if(j <= n && j > i){
#      temp <- matrix[j,i] / pivot
#      matrix[j,i] <- temp
#      matrix[j,j] = matrix[j,j] - temp*temp*pivot
#      for(k in seq(j+1, n)){
#       if(k <= n && k > j){
#        matrix[k, j] <- matrix[k, j] - temp * matrix[k, i]
#       }
#      }
#     }
#    }
#   }
#  }
#  list(matrix = matrix,
#       rank = rank * nonneg)
# }
#
# chsolve2 <- function(matrix, y){
#
#  n <- ncol(matrix)
#
#  for (i in seq(n)) {
#   temp = y[i]
#   for (j in seq(i)){
#    if(j < i)
#     temp = temp - y[j] * matrix[i, j]
#   }
#   y[i] = temp
#  }
#
#  for (i in seq(n,1)){
#   if (matrix[i, i]==0){
#    y[i] =0;
#   } else {
#    temp = y[i] / matrix[i, i]
#    for (j in seq(i+1, n)){
#     if(j <= n && j > i){
#      temp = temp - y[j] * matrix[j, i]
#     }
#    }
#    y[i] = temp
#   }
#  }
#  y
# }
#
# chinv2 <- function(matrix){
#
#  n <- ncol(matrix)
#
#  for (i in seq(n)){
#
#   if (matrix[i,i] >0){
#    matrix[i,i] = 1/matrix[i,i]
#
#    for (j in seq(i+1, n)) {
#     if(j <= n && j > i){
#      matrix[j,i] = -matrix[j,i]
#      for (k in seq(i)){
#       if(k < i)
#        matrix[j,k] = matrix[j,k] + matrix[j,i]*matrix[i,k]
#      }
#     }
#    }
#   }
#  }
#
#  for (i in seq(n)) {
#   if (matrix[i,i]==0) {
#    for (j in seq(i)) if (j < i) matrix[i,j]=0;
#    for (j in seq(i, n)) if (j < i) matrix[j,i]=0;
#   } else {
#    for (j in seq(i+1, n)) {
#     if (j > i && j <= n){
#      temp = matrix[j,i] * matrix[j,j];
#      if (j!=i) matrix[i,j] = temp;
#      for (k in seq(i,j)){
#       if (k < j){
#        matrix[i,k] = matrix[i,k] + temp*matrix[j,k];
#       }
#      }
#     }
#    }
#   }
#  }
#
#  matrix
#
# }
#
# cph_iter <- function(X,
#                      beta,
#                      y_time,
#                      y_status,
#                      weights,
#                      means,
#                      scales){
#
#  # u contains score vector
#  u <- rep(0, nvar)
#  # a contains weighted sums of x
#  a <- rep(0, nvar)
#  # a2 reserved for observations with events
#  a2 <- rep(0, nvar)
#  # what is imat?
#  imat <- matrix(0, ncol = nvar, nrow = nvar)
#  # cmat contains sums of squares
#  cmat <- matrix(0, ncol = nvar, nrow = nvar)
#  # cmat2 reserved for observations with events
#  cmat2 <- matrix(0, ncol = nvar, nrow = nvar)
#
#  nrisk  = 0 # number at risk for death at this time point
#  denom  = 0 # sum of risks for the deaths
#  loglik = 0 # log likelihood
#
#  # go from last person's time to first person's time
#  # to make the cumulative sum easier to manage
#  person = nrow(X)
#
#  while(person >= 1){
#
#   person_time = y_time[person] # time of event for current person
#   ndead  = 0 # number of deaths at this time point
#   deadwt = 0 # sum of weights for the deaths
#   denom2 = 0 # sum of weighted risks for the deaths
#
#   while(person >= 1 && y_time[person] == person_time){
#
#    nrisk <- nrisk + 1
#
#    # form the dot product of X[person,] and beta
#    x_beta <- 0
#
#    for(i in seq(nvar)){
#     x_beta <- x_beta + X[person, i] * beta[i]
#    }
#
#    risk <- exp(x_beta) * weights[person]
#
#    if(y_status[person] == 0){
#
#     denom = denom + risk
#
#     for (i in seq(nvar)) {
#
#      a[i] = a[i] + risk * X[person, i]
#
#      for (j in seq(i)){
#       cmat[j,i] = cmat[j,i] + risk * X[person, i] * X[person, j]
#      }
#
#     }
#
#    }
#
#    if (y_status[person] == 1){
#
#     ndead  <- ndead  + 1
#     deadwt <- deadwt + weights[person]
#     denom2 <- denom2 + risk
#     loglik <- loglik + weights[person] * x_beta
#
#     for (i in seq(nvar)) {
#
#      u[i]  <-  u[i] + weights[person] * X[person, i];
#      a2[i] <-  a2[i] + risk * X[person, i]
#
#      for (j in seq(i)){
#       cmat2[j,i] <- cmat2[j,i] + risk * X[person, i] * X[person, j]
#      }
#
#     }
#
#    }
#
#    person <- person - 1
#
#   }
#
#   # if number of events at this time is >0,
#   # then we add to the main terms:
#   if(ndead > 0){
#
#    # breslow method for ties
#    if(method == 0 || ndead == 1){
#
#     denom <- denom + denom2;
#     loglik <- loglik - deadwt * log(denom);
#
#     for (i in seq(nvar)) {
#      a[i]  = a[i] + a2[i]
#      temp2 = a[i] / denom  # Mean
#      u[i]  = u[i] - deadwt * temp2
#
#      for (j in seq(i)) {
#       cmat[j,i] = cmat[j,i] + cmat2[j,i]
#       imat[j,i] = imat[j,i] + deadwt * (cmat[j,i] - temp2 * a[j]) / denom;
#      }
#     }
#
#    } else {
#
#     # Efron
#     #
#     # If there are 3 deaths we have 3 terms: in the first the
#     #  three deaths are all in, in the second they are 2/3
#     #  in the sums, and in the last 1/3 in the sum.  Let k go
#     #  1 to ndead: we sequentially add a2/ndead and cmat2/ndead
#     #  and efron_wt/ndead to the totals.
#     #
#     wtave = deadwt/ndead;
#
#     for (k in seq(ndead)) {
#
#      denom = denom + denom2/ndead;
#      loglik = loglik - wtave * log(denom);
#
#      for (i in seq(nvar)) {
#
#       a[i]  = a[i] + a2[i] / ndead;
#       temp2 = a[i] / denom;
#       u[i]  = u[i] - wtave * temp2;
#
#       for (j in seq(i)) {
#
#        cmat[j,i] = cmat[j,i] + cmat2[j,i]/ndead;
#        imat[j,i] = imat[j,i] + wtave*(cmat[j,i] - temp2*a[j]) / denom;
#
#       }
#
#      }
#
#     }
#
#    }
#
#    for (i in seq(nvar)) {
#     a2[i] = 0
#     for (j in seq(nvar)) cmat2[j, i] = 0
#    }
#
#   }
#
#  }
#
#  list(imat = imat,
#       loglik = loglik,
#       u = u)
#
# }
