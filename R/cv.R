#' Cross-validation for multiLMMsel
#'
#' Perform cross-validation over a sequence of \code{lambda}
#' values (per response) and \code{tau} values to select tuning parameters for
#' the multiLMMsel algorithm.
#'
#' The function first selects the random-effects structure for each response,
#' then refits the covariance matrices and finally tunes the fixed-effects
#' penalty via cross-validation.
#'
#' @param Y Response matrix of size \code{N x d}.
#' @param Z Random-effects design matrix of size \code{N x q}.
#' @param X Fixed-effects design matrix of size \code{N x p}.
#' @param lambda.path Numeric vector of candidate regularization parameters
#'   for random effects selection.
#' @param tau.path Numeric vector of candidate regularization parameters
#'   for fixed effects selection.
#' @param threshold Threshold on the diagonal of \code{Sigma_B} for defining
#'   the active set.
#' @param folds Number of CV folds at the subject/cluster level.
#' @param id Subject/cluster index for each row (same length as \code{nrow(Y)}).
#' @param gamma Noise variance parameter used in the bootstrap-based beta
#'   updating scheme.
#' @param gamma.weight Power used in the adaptive weights for the diagonal of
#'   \code{Sigma_B}.
#' @param alpha One-standard-error parameter deciding how much standard
#'   error to use to pick a more parsimonious lambda along the CV curve.
#' @param sigmaB,sigmaE,beta Optional "true" values used when evaluating F-norm
#'   and F1 scores in simulation studies. If all are \code{NULL}, only estimates
#'   are returned.
#' @param eta Step size for the update of \code{Sigma_B}.
#'   If \code{NULL}, it is chosen theoretically from the largest eigenvalue of the Hessian.
#' @param max.iterB,tolB Parameters for the inner SigmaB update.
#' @param bootstrap.iter Number of bootstrap iterations used for beta estimation.
#'
#' @return If \code{sigmaB}, \code{sigmaE}, and \code{beta} are \code{NULL},
#'   a list with:
#'   \item{lambda.optimal}{Vector of selected lambda for each response.}
#'   \item{tau.optimal}{Selected tau.}
#'   \item{l.lambda}{Matrix of CV likelihoods over lambda.}
#'   \item{l.tau}{Matrix of CV likelihoods over tau.}
#'   \item{SigmaB}{Estimated random-effects covariance.}
#'   \item{SigmaE}{Estimated error covariance.}
#'   \item{beta}{Estimated fixed-effects matrix.}
#'
#'   Otherwise, additional F-norm and F1-score summaries are returned for
#'   simulation evaluation.
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#' d <- 3
#' q <- 5
#' dq <- d*q
#' p <- 10
#' m <- 50
#' ni <- 6
#' s <- 3
#' rho <- 0.5
#' M <- bdiag(toeplitz(rho^seq(0, s-1)),matrix(0,q-s,q-s))
#' R <- cbind(M,M,M)
#' sigmaB <- rbind(R,R,R)
#' rho <- 0.75
#' sigmaE <- toeplitz(rho^seq(0, d - 1))
#' id <- rep(1:m,each=ni)
#' m <- length(unique(id)) # Number of subjects (groups)
#' N <- length(id) # Total number of observations
#' epsilon <- mvrnorm(N,rep(0,d),sigmaE)
#' vecB <- mvrnorm(m,rep(0,dq),sigmaB)
#' B <- do.call(rbind,lapply(1:nrow(vecB), function(i) {
#'   matrix(vecB[i, ], nrow = q, ncol = d, byrow = FALSE)
#' }))
#' Z <- matrix(rnorm(N*q),nrow = N,ncol = q)
#' X <- matrix(rnorm(N*p),nrow = N,ncol = p)
#' beta <- matrix(c(1,-1,2,rep(0,7),-1,2,-1,rep(0,7),1,1,2,rep(0,7)),nrow = p,ncol = d)
#' Z_list <- lapply(split(seq_len(nrow(Z)), id),
#'                  function(idx) Z[idx, , drop = FALSE])
#' Z_long <- bdiag(lapply(Z_list, as.matrix))
#' Y <- X %*% beta + as.matrix(Z_long) %*% B + epsilon
#'
#' lambda.max <- 5000
#' lambda.min <- lambda.max/100
#' length.out <- 30
#' lambdas <- exp(seq(log(lambda.max),log(lambda.min),length.out=length.out))
#'
#' tau.max <- 1
#' tau.min <- tau.max/100
#' length.out <- 30
#' taus <- seq(tau.max,tau.min,length.out=length.out)
#'
#' result.sim <- multiLMMsel.cv(
#'   Y,Z,X,lambdas,taus,id = id,
#'   sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,alpha = 1
#'   )
#' }
#'
#' @useDynLib multiLMMsel, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import Matrix
#' @importFrom RSpectra eigs
#' @importFrom glmnet glmnet
#' @importFrom stats rnorm sd
#' @importFrom methods as
#' @export
multiLMMsel.cv <- function(Y,Z,X,lambda.path,tau.path,threshold=0.01,folds=5,id,gamma=2,gamma.weight=2,alpha=1,sigmaB=NULL,sigmaE=NULL,beta=NULL,eta=NULL,max.iterB=1000,tolB=1e-4,bootstrap.iter=30){
  id <- as.character(id)
  ord <- order(match(id, unique(id)))
  id <- id[ord]
  Y <- as.matrix(Y)[ord, , drop = FALSE]
  X <- as.matrix(X)[ord, , drop = FALSE]
  Z <- as.matrix(Z)[ord, , drop = FALSE]
  id <- factor(id, levels = unique(id))
  nis <- tabulate(id)
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z <- bdiag(lapply(Z_list, as.matrix))

  N <- sum(nis)
  m <- length(nis)
  d <- ncol(Y)
  q <- ncol(Z)/m
  p <- ncol(X)
  test.size <- floor(m/folds)
  lambda.size <- length(lambda.path)
  tau.size <- length(tau.path)
  l.lambda <- matrix(nrow = lambda.size, ncol = folds)
  l.tau <- matrix(nrow = tau.size, ncol = folds)

  beta.d <- solve(t(X) %*% X, t(X) %*% Y)
  y.trisum.d <- triple_sum_y(as.matrix(Y),as.matrix(Z),X,as.matrix(beta.d),nis,m,d,q)
  ww.trisum.d <- triple_sum_ww(as.matrix(Z),nis,m,d,q)
  diag.vec <- diag(updateSigmaB(sigmaB=NULL,eta=NULL,ww.trisum.d,y.trisum.d,lambda.i=0,Y,Z,d,q,delta=0,max.iterB,tolB)$sigmaB)

  lambda.lmax.list <- numeric(d)
  lambda.optimal.list <- numeric(d)
  ww.trisum <- NULL
  eta <- NULL
  active.set.vec <- NULL
  active.set.mat <- NULL
  ww.trisum.all <- triple_sum_ww(as.matrix(Z),nis,m,d=1,q)

  for (l in 1:d){
    Y_l <- as.matrix(Y[,l])
    diag_l <- diag.vec[(1+(l-1)*q):(l*q)]
    weight_l <- 1/(diag_l^gamma.weight)
    for (fold in 1:folds) {
      if (fold == 1){
        y.test <- Y_l[1:cumsum(nis)[test.size],]
        y.train <- Y_l[-c(1:cumsum(nis)[test.size]),]
        x.test <- X[1:cumsum(nis)[test.size],]
        x.train <- X[-c(1:cumsum(nis)[test.size]),]
        z.test <- Z[1:cumsum(nis)[test.size], ((fold-1)*test.size*q+1):(fold*test.size*q)]
        z.train <- Z[-c(1:cumsum(nis)[test.size]), -c(((fold-1)*test.size*q+1):(fold*test.size*q))]
      }
      else {
        y.test <- Y_l[(cumsum(nis)[(fold-1)* test.size]+1):cumsum(nis)[fold*test.size],]
        y.train <- Y_l[-c((cumsum(nis)[(fold-1)* test.size]+1):cumsum(nis)[fold*test.size]),]
        x.test <- X[(cumsum(nis)[(fold-1)* test.size]+1):cumsum(nis)[fold*test.size],]
        x.train <- X[-c((cumsum(nis)[(fold-1)* test.size]+1):cumsum(nis)[fold*test.size]),]
        z.test <- Z[(cumsum(nis)[(fold-1)* test.size]+1):cumsum(nis)[fold*test.size], ((fold-1)*test.size*q+1):(fold*test.size*q)]
        z.train <- Z[-c((cumsum(nis)[(fold-1)* test.size]+1):cumsum(nis)[fold*test.size]), -c(((fold-1)*test.size*q+1):(fold*test.size*q))]
      }

      nis.test <- nis[((fold-1)*test.size+1):(fold*test.size)]
      nis.train <- nis[-c(((fold-1)*test.size+1):(fold*test.size))]
      m.train <- length(nis.train)
      bbeta <- solve(t(as.matrix(x.train)) %*% as.matrix(x.train), t(as.matrix(x.train)) %*% as.matrix(y.train))
      if (l == 1){
        ww.trisum[[fold]] <- triple_sum_ww(as.matrix(z.train),nis.train,m.train,d=1,q)
        if (nrow(ww.trisum[[fold]]) <= 2){
          L <- norm(ww.trisum[[fold]],"2")
        }
        else{
          L <- eigs(ww.trisum[[fold]], k = 1, which = "LM")$values
        }
        eta[[fold]] <- 1/L
      }
      y.trisum <- triple_sum_y(as.matrix(y.train),as.matrix(z.train),as.matrix(x.train),bbeta,nis.train,m.train,d=1,q)
      y.dbsum <- double_sum_y(as.matrix(y.train),as.matrix(x.train),bbeta,nis.train,m.train,d=1)

      sol.path <- sigma.path(as.matrix(y.train),as.matrix(z.train),as.matrix(x.train),lambda.path,eta[[fold]],nis.train,d=1,q,
                             weight=weight_l,ww.trisum[[fold]],y.trisum,y.dbsum,max.iterB=max.iterB,tolB=tolB)
      sigmaB.path <- sol.path$sigmaB.path
      sigmaE.path <- sol.path$sigmaE.path

      for (lambda.index in 1:length(lambda.path)) {
        l.lambda[lambda.index,fold] <- likelihood(as.matrix(y.test),as.matrix(z.test),as.matrix(x.test),
                                                  sigmaB.path[[lambda.index]],sigmaE.path[[lambda.index]],bbeta,nis.test,d=1,q)
      }
      cat("Fold", fold, "lambda completed\n")
    }

    l.mean <- apply(l.lambda, 1, mean)
    l.max <- max(l.mean)
    l.max.index <- which.max(l.mean)
    lambda.lmax <- lambda.path[l.max.index]
    lambda.lmax.list[l] <- lambda.lmax

    l.se <- sd(l.lambda[l.max.index,])/sqrt(folds)
    l.optimal.index <- which(l.mean >= l.max - alpha * l.se)[1]
    lambda.optimal <- lambda.path[l.optimal.index]
    lambda.optimal.list[l] <- lambda.optimal

    lambda.weight.l <- lambda.optimal * diag(weight_l)
    y.trisum.all <- triple_sum_y(as.matrix(Y_l),as.matrix(Z),as.matrix(X),as.matrix(beta.d[,l]),nis,m,d=1,q)
    sigmaB.select <- updateSigmaB(sigmaB.path[[l.optimal.index]],eta=NULL,ww.trisum.all,y.trisum.all,lambda.weight.l,Y_l,Z,d=1,q,delta=0,max.iterB,tolB)
    sigmaB_l <- sigmaB.select$sigmaB
    active.set.mat[[l]] <- which(diag(sigmaB_l) > threshold)
    active.set.vec <- c(active.set.vec,(active.set.mat[[l]]+q*(l-1)))
  }

  #Estimate Sigma_B
  n <- length(active.set.vec)
  sigmaB.hat <- estimateSigmaB(sigmaB=NULL,eta=NULL,Y,X,Z,beta.d,delta=0,active.set.vec,active.set.mat,n,nis,m,d,q)$sigmaB

  #Estimate Sigma_e
  y.dbsum.all <- double_sum_y(Y,X,beta.d,nis,m,d)
  sigmaE.hat <- updateSigmaE(y.dbsum.all,Z,sigmaB.hat,nis,m,d,N,q,delta=1e-4)

  #cv for beta
  for (fold in 1:folds) {
    if (fold == 1){
      y.test <- Y[1:cumsum(nis)[test.size],]
      y.train <- Y[-c(1:cumsum(nis)[test.size]),]
      x.test <- X[1:cumsum(nis)[test.size],]
      x.train <- X[-c(1:cumsum(nis)[test.size]),]
      z.test <- Z[1:cumsum(nis)[test.size], ((fold-1)*test.size*q+1):(fold*test.size*q)]
      z.train <- Z[-c(1:cumsum(nis)[test.size]), -c(((fold-1)*test.size*q+1):(fold*test.size*q))]
    }
    else {
      y.test <- Y[(cumsum(nis)[(fold-1)* test.size]+1):cumsum(nis)[fold*test.size],]
      y.train <- Y[-c((cumsum(nis)[(fold-1)* test.size]+1):cumsum(nis)[fold*test.size]),]
      x.test <- X[(cumsum(nis)[(fold-1)* test.size]+1):cumsum(nis)[fold*test.size],]
      x.train <- X[-c((cumsum(nis)[(fold-1)* test.size]+1):cumsum(nis)[fold*test.size]),]
      z.test <- Z[(cumsum(nis)[(fold-1)* test.size]+1):cumsum(nis)[fold*test.size], ((fold-1)*test.size*q+1):(fold*test.size*q)]
      z.train <- Z[-c((cumsum(nis)[(fold-1)* test.size]+1):cumsum(nis)[fold*test.size]), -c(((fold-1)*test.size*q+1):(fold*test.size*q))]
    }
    nis.test <- nis[((fold-1)*test.size+1):(fold*test.size)]
    nis.train <- nis[-c(((fold-1)*test.size+1):(fold*test.size))]
    bbeta <- solve(t(as.matrix(x.train)) %*% as.matrix(x.train), t(as.matrix(x.train)) %*% as.matrix(y.train))

    bbeta.path <- beta.path(as.matrix(y.train),as.matrix(z.train),as.matrix(x.train),tau.path,gamma,sigmaB.hat,sigmaE.hat,bbeta,nis.train,d,p,q,bootstrap.iter)

    for (tau.index in 1:length(tau.path)) {
      l.tau[tau.index,fold] <- likelihood(as.matrix(y.test),as.matrix(z.test),as.matrix(x.test),
                                          sigmaB.hat,sigmaE.hat,bbeta.path[[tau.index]],nis.test,d,q)
    }
    cat("Fold", fold, "tau completed\n")
  }

  l.max.index <- which.max(apply(l.tau, 1, mean))
  tau.optimal <- tau.path[l.max.index]

  beta.hat <- updateBeta(Y,X,Z,tau.optimal,gamma,sigmaB.hat,sigmaE.hat,nis,d,p,q,N,bootstrap.iter)

  if ((is.null(sigmaB) == TRUE) | (is.null(sigmaE) == TRUE) | (is.null(beta) == TRUE)){
    return(list(lambda.optimal=lambda.optimal.list,tau.optimal=tau.optimal,l.lambda=l.lambda,l.tau=l.tau,SigmaB=sigmaB.hat,SigmaE=sigmaE.hat,beta=beta.hat))
  }
  else{
    fnorm.B <- norm(sigmaB-sigmaB.hat,"F")
    #onorm.B <- norm(sigmaB-sol$SigmaB,"2")

    diagB.true <- diag(sigmaB) != 0
    diagB.pred  <-  diag(sigmaB.hat) != 0
    #TPR <- sum(diagB.pred == 1 & diagB.true == 1) / sum(diagB.true==1)
    #FPR <- sum(diagB.pred == 1 & diagB.true == 0) / sum(diagB.true==0)
    #TNR <- sum(diagB.pred == 0 & diagB.true == 0) / sum(diagB.true==0)
    #FNR <- sum(diagB.pred == 0 & diagB.true == 1) / sum(diagB.true==1)
    #FDP <- sum(diagB.pred == 1 & diagB.true == 0) / max(sum(diagB.pred==1),1)
    F1.B <- 2*sum(diagB.pred == 1 & diagB.true == 1)/(2*sum(diagB.pred == 1 & diagB.true == 1)+sum(diagB.pred == 1 & diagB.true == 0)+sum(diagB.pred == 0 & diagB.true == 1))

    fnorm.beta <- norm(beta-beta.hat,"F")
    beta.true <- as.vector(beta) != 0
    beta.pred <- as.vector(beta.hat) != 0
    F1.beta <- 2*sum(beta.pred == 1 & beta.true == 1)/(2*sum(beta.pred == 1 & beta.true == 1)+sum(beta.pred == 1 & beta.true == 0)+sum(beta.pred == 0 & beta.true == 1))


    return(list(fnorm.B=fnorm.B,F1.B=F1.B,fnorm.beta=fnorm.beta,F1.beta=F1.beta,
                l.lambda=l.lambda,l.tau=l.tau,lambda.optimal=lambda.optimal.list,tau.optimal=tau.optimal,
                SigmaB=sigmaB.hat,SigmaE=sigmaE.hat,beta=beta.hat,true.sigmaB=sigmaB,true.sigmaE=sigmaE,true.beta=beta))
  }
}
