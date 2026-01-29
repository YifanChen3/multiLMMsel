# Objective function when minimizing Sigma_B
objB <- function(Y,X,Z,sigmaB,lambda,nis,m,d,q){
  result <- NULL
  for (i in 1:m) {
    for (j in 1:nis[i]) {
      for (k in 1:nis[i]) {
        if (j != k){
          each <- norm((kronecker(Diagonal(d),t(Z[zij(i,j,nis,q)[[1]],zij(i,j,nis,q)[[2]]])) %*% sigmaB %*% kronecker(Diagonal(d),Z[zij(i,k,nis,q)[[1]],zij(i,k,nis,q)[[2]]]) - (Y[yij(i,j,nis),]-t(beta) %*% X[yij(i,j,nis),]) %*% t(Y[yij(i,k,nis),]-t(beta)%*%X[yij(i,k,nis),])),"F")
          result <- c(result, each^2)
        }
      }
    }
  }
  return(result)
}


# Project A onto a Positive Semi-definite Cone
prox.psd <- function(A,delta){
  eigenA <- eigen(A)
  U <- eigenA$vectors
  if (length(U) == 1){
    L <- pmax(eigenA$values,delta)
    return(as.matrix(L))
  }
  else{
    L <- diag(pmax(eigenA$values,delta))
    result <- U %*% L %*% t(U)
    return(result)
  }
}

prox.l1 <- function(A,eta,lambda){
  diag(A) <- pmax(diag(A)-eta*lambda,0)
  return(A)
}

# Find (i,j)-th index of y_ij from the big matrix Y
yij <- function(i,j,nis){
  if (i == 1) {
    result <- j
  }
  else {
    result <- cumsum(nis)[i-1] + j
  }
  return(result)
}

# Find (i,j)-th index of z_ij from the big matrix Z
zij <- function(i,j,nis,q){
  if (i == 1) {
    result.i <- j
  }
  else {
    result.i <- cumsum(nis)[i-1] + j
  }
  start.j <- q*(i-1) + 1
  end.j <- q*i
  result.j <- c(start.j:end.j)
  return(list(result.i,result.j))
}

# Combine together
gradient.loss <- function(ww.trisum,sigmaB,y.trisum,lambda.i,d,q){
  result <- matrix(ww.trisum %*% as.vector(sigmaB),d*q,d*q) + y.trisum + lambda.i
  return(result)
}

# Gradient Descent for 1 step
sigmaB.tilde <- function(sigmaB,eta,ww.trisum,y.trisum,lambda.i,d,q){
  result <- sigmaB - eta * gradient.loss(ww.trisum,sigmaB,y.trisum,lambda.i,d,q)
  return(result)
}

updateSigmaE <- function(y.dbsum,Z,sigmaB,nis,m,d,N,q,delta){
  result <- prox.psd((y.dbsum-double_sum_z_b(as.matrix(Z),as.matrix(sigmaB),nis,m,d,q))/N,delta)
  return(result)
}

updateSigmaB <- function(sigmaB=NULL,eta=NULL,ww.trisum,y.trisum,lambda.i,Y,Z,d,q,delta,max.iter,tol){
  iter <- 0
  if (is.null(sigmaB) == TRUE){
    sigmaB <- Diagonal(d*q)
  }
  if (is.null(eta) == TRUE){
    if (nrow(ww.trisum) <= 2){
      L <- norm(ww.trisum,"2")
    }
    else{
      L <- eigs(ww.trisum, k = 1, which = "LM")$values
    }
    eta <- 1/L
  }

  residual.list <- NULL
    while (iter < max.iter) {
      iter <- iter + 1
      sigmaB.old <- sigmaB
      sigmaB.pre.prox <- sigmaB.tilde(sigmaB,eta,ww.trisum,y.trisum,lambda.i,d,q)
      sigmaB <- prox.psd(sigmaB.pre.prox,delta)
      residual <- norm(sigmaB-sigmaB.old,"F")
      residual.list <- c(residual.list,residual)
      if (residual/max(1,residual) < tol){
        break
      }
    }
    return(list(sigmaB=sigmaB,iter=iter,residual.list=residual.list))
}

gradient.est <- function(sigmaB,z4.trisum,zzy.trisum,n){
  result <- matrix(z4.trisum %*% as.vector(sigmaB),n,n) + zzy.trisum
  return(result)
}

sigmaB.est.gd <- function(sigmaB,eta,z4.trisum,zzy.trisum,n){
  result <- sigmaB - eta * gradient.est(sigmaB,z4.trisum,zzy.trisum,n)
  return(result)
}


estimateSigmaB <- function(sigmaB=NULL,eta=NULL,Y,X,Z,beta,delta,active.set.vec,active.set.mat,n,nis,m,d,q,max.iter=1000,tol=1e-8){
  iter <- 0
  if (n == 0){
    stop("The tunning parameter lambda is too large, so the model is penalizing all the random effects to zero. Please use a smaller lambda.")
  }
  if (is.null(sigmaB) == TRUE){
    sigmaB <- Diagonal(n)
  }
  active.set.mat.c <- lapply(c(1:d), function(i){active.set.mat[[i]]-1}) #Change to 0-base for Cpp
  two.trisums <- triple_sum_for_estimate(Y,X,as.matrix(Z),beta,nis,m,d,q,active.set.mat.c,n)
  z4.trisum <- two.trisums$z4
  zzy.trisum <- two.trisums$zzy

  if (is.null(eta) == TRUE){
    if (nrow(z4.trisum) <= 2){
      L <- norm(z4.trisum,"2")
    }
    else{
      L <- eigs(z4.trisum, k = 1, which = "LM")$values
    }
    eta <- 1/L
  }
  while (iter < max.iter) {
    iter <- iter + 1
    sigmaB.old <- sigmaB
    sigmaB.pre.prox <- sigmaB.est.gd(sigmaB,eta,z4.trisum,zzy.trisum,n)
    sigmaB <- prox.psd(sigmaB.pre.prox,delta)
    residual <- norm(sigmaB-sigmaB.old,"F")
    if (residual/max(1,residual) < tol){
      break
    }
  }

  #impute 0's
  full.matrix <- matrix(0,d*q,d*q)
  full.matrix[active.set.vec,active.set.vec] <- sigmaB
  return(list(sigmaB=full.matrix,iter=iter))
}


vec.yi.t <- function(Y,nis,i){
  if (i == 1){
    index <- 1
  }
  else {
    index <- cumsum(nis)[i-1] + 1
  }
  Yi <- Y[index:(index+nis[i]-1),]
  result <- as.vector(t(Yi))
  return(result)
}

sigmaT <- function(Z,sigmaB,sigmaE,nis,d,q){
  m <- length(nis)
  sigma.list <- vector("list",m)
  for (i in 1:m) {
    sigma.list[[i]] <- cov_yi(as.matrix(Z),as.matrix(sigmaB),sigmaE,nis,i,d,q)
  }
  result <- as.matrix(bdiag(sigma.list))
  return(result)
}

updateBeta <- function(Y,X,Z,tau,gamma,sigmaB,sigmaE,nis,d,p,q,N,bootstrap.iter=30){
  XI <- kronecker(X,Diagonal(d))
  sigma.t.base <- sigmaT(Z,sigmaB,sigmaE,nis,d,q)
  sigma.t.base <- as(sigma.t.base,"dgCMatrix")
  result <- matrix(0,nrow = d,ncol = p)
  select.count <- numeric(d * p)

  for (i in 1:bootstrap.iter) {
    sigma.t <- sigma.t.base + gamma * diag(N*d) #Updated here
    vec.Y.t <- as.vector(t(Y)) + rnorm(N*d,0,sqrt(gamma)) #Updated here, also the parameter of this function
    U <- chol(sigma.t)
    x.tilda <- solve(t(U), XI)
    y.tilda <- solve(t(U), vec.Y.t)
    beta.vec <- solve(crossprod(x.tilda), crossprod(x.tilda, y.tilda))

    weights <- 1/(abs(beta.vec)+1e-8)
    model <- glmnet(x.tilda,y.tilda,alpha = 1,lambda = tau,penalty.factor = weights)
    beta.hat <- model$beta
    select.count[which(beta.hat != 0)] <- select.count[which(beta.hat != 0)] + 1
    result <- result + matrix(beta.hat, nrow = d)
  }
  result <- result/bootstrap.iter
  support <- which(select.count >= (bootstrap.iter - 1))
  result[setdiff(c(1:(d*p)),support)] <- 0

  return(t(result))
}

likelihood <- function(Y,Z,X,sigmaB,sigmaE,beta,nis,d,q){
  N <- sum(nis)
  m <- length(nis)
  result <- N * d * log(2 * pi)

  sigma.t <- sigmaT(Z,sigmaB,sigmaE,nis,d,q)
  sigma.t <- as(sigma.t,"dgCMatrix")
  result <- result + determinant(sigma.t, logarithm = TRUE)$modulus[1]

  y.mu <- matrix(as.vector(t(Y)) - kronecker(X,Diagonal(d)) %*% as.vector(t(beta)),ncol=1)
  pdf.kernel <- t(y.mu) %*% solve(sigma.t,y.mu)
  result <- result + pdf.kernel

  result <- -0.5 * result

  return(as.numeric(result))
}





sigma.path <- function(Y,Z,X,lambda.path,eta,nis,d,q,weight,ww.trisum,y.trisum,y.dbsum,max.iterB=1000,tolB=1e-8){
  N <- sum(nis)
  m <- length(nis)
  dq <- d*q

  sigmaB.path <- NULL
  sigmaE.path <- NULL
  sigmaB.path[[1]] <- matrix(0,dq,dq)
  sigmaE.path[[1]] <- prox.psd(y.dbsum/N,1e-4)
  iterB.path <- NULL
  index <- 0
  for (lambda in lambda.path) {
    index <- index + 1
    lambda.weight <- lambda * diag(weight)
    if (index == 1){
      sigmaB.sol <- updateSigmaB(sigmaB.path[[1]],eta,ww.trisum,y.trisum,lambda.weight,Y,Z,d=1,q,delta=0,max.iterB,tolB)
    }
    else{
      sigmaB.sol <- updateSigmaB(sigmaB.path[[index-1]],eta,ww.trisum,y.trisum,lambda.weight,Y,Z,d=1,q,delta=0,max.iterB,tolB)
    }
    sigmaB <- sigmaB.sol$sigmaB
    sigmaE <- updateSigmaE(y.dbsum,Z,sigmaB,nis,m,d,N,q,delta=1e-4)
    iterB.path <- c(iterB.path,sigmaB.sol$iterB)
    sigmaB.path[[index]] <- sigmaB
    sigmaE.path[[index]] <- sigmaE
  }
  return(list(sigmaB.path=sigmaB.path,sigmaE.path=sigmaE.path,iterB.path=iterB.path))
}

beta.path <- function(Y,Z,X,tau.path,gamma,sigmaB,sigmaE,beta,nis,d,p,q,bootstrap.iter=30){
  N <- sum(nis)
  bbeta.path <- NULL
  bbeta.path[[1]] <- beta

  XI <- kronecker(X,Diagonal(d))
  sigma.t.base <- sigmaT(Z,sigmaB,sigmaE,nis,d,q)
  sigma.t.base <- as(sigma.t.base,"dgCMatrix")
  beta.matrix <- matrix(0,nrow = d*p,ncol = length(tau.path))

  for (i in 1:bootstrap.iter) {
    sigma.t <- sigma.t.base + gamma * diag(N*d)
    vec.Y.t <- as.vector(t(Y)) + rnorm(N*d,0,sqrt(gamma))
    U <- chol(sigma.t)
    x.tilda <- solve(t(U), XI)
    y.tilda <- solve(t(U), vec.Y.t)
    beta.vec <- solve(crossprod(x.tilda), crossprod(x.tilda, y.tilda))

    weights <- 1/(abs(beta.vec)+1e-8)
    model <- glmnet(x.tilda,y.tilda,alpha = 1,lambda = tau.path,penalty.factor = weights)
    beta.matrix <- beta.matrix + model$beta
  }
  beta.matrix <- beta.matrix/bootstrap.iter


  index <- 0
  for (tau in tau.path) {
    index <- index + 1
    bbeta <- t(matrix(beta.matrix[,index],nrow = d))
    bbeta.path[[index]] <- bbeta
  }
  return(bbeta.path)
}


