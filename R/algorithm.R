#' Variable Selection for Multivariate Mixed-Effects Models (multiLMMsel)
#'
#' Fit a multivariate high-dimensional mixed-effects model with variable
#' selection using the proposed moment-based convex optimization algorithm.
#'
#' The function performs five-step updates of the random-effects covariance
#' \eqn{\Sigma_B}, the residual covariance \eqn{\Sigma_\varepsilon}, and the
#' fixed effects \eqn{\beta}, given user-specified regularization parameters
#' \code{lambdas} and \code{tau}.
#'
#' @param Y Numeric matrix of responses of size \code{N x d}.
#' @param Z Random-effects design matrix of size \code{N x q}.
#' @param X Fixed-effects design matrix of size \code{N x p}.
#' @param lambdas Numeric vector of length \code{d} giving the penalty level
#'   for each response dimension.
#' @param tau Regularization parameter for the adaptive lasso step in the
#'   fixed-effects update.
#' @param threshold Threshold on the diagonal of \code{Sigma_B} used to define
#'   the active set of random effects.
#' @param gamma Noise variance parameter used in the bootstrap-based beta
#'   updating scheme.
#' @param gamma.weight Power used in the adaptive weights for the diagonal of
#'   \code{Sigma_B}.
#' @param eta Step size for the proximal-gradient update of \code{Sigma_B}.
#'   If \code{NULL}, it is chosen theoretically from the largest eigenvalue of the Hessian.
#' @param d Number of responses. If \code{NULL}, taken as \code{ncol(Y)}.
#' @param p Number of fixed effects. If \code{NULL}, taken as \code{ncol(X)}.
#' @param q Number of random effects per response. If \code{NULL}, inferred
#'   from \code{ncol(Z)}.
#' @param id Integer or factor vector of length \code{N} giving the subject or
#'   cluster index for each row of \code{Y}, \code{X}, and \code{Z}.
#' @param deltaB Eigenvalue threshold for the projection of \code{Sigma_B} onto PSD cone.
#' @param deltaE Eigenvalue threshold for the projection of \code{Sigma_E} onto PSD cone.
#' @param max.iter Maximum number of outer iterations.
#' @param tol Convergence tolerance outer iterations.
#' @param max.iterB Maximum number of inner iterations for \code{Sigma_B}
#'   updates.
#' @param tolB Convergence tolerance for inner \code{Sigma_B} updates.
#' @param sigmaB Optional initial value for \code{Sigma_B}.
#' @param beta Optional initial value for the fixed effects matrix \code{beta}.
#' @param y.dbsum,ww.trisum,y.trisum,ww.trisum.d,y.trisum.d Optional
#'   precomputed quantities (typically produced by C++ helper
#'   functions) to avoid recomputation.
#' @param bootstrap.iter Number of bootstrap iterations used in the fixed-effect
#'   update.
#'
#' @return A list with components:
#'   \item{SigmaB}{Estimated random-effects covariance matrix.}
#'   \item{SigmaE}{Estimated residual covariance matrix.}
#'   \item{beta}{Estimated fixed-effects coefficient matrix.}
#'   \item{iterB}{Number of iterations for SigmaB updates.}
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
#' lambdas <- c(200,200,200)
#' tau <- 0.2
#' result <- multiLMMsel(Y,Z,X,lambdas,tau,id = id)
#' }
#'
#' @export
multiLMMsel <- function(Y,Z,X,lambdas,tau,threshold=0.01,gamma=2,gamma.weight=2,eta=NULL,d=NULL,p=NULL,q=NULL,id,deltaB=0,deltaE=1e-4,
                  max.iter=1,tol=1e-2,max.iterB=1000,tolB=1e-8,
                  sigmaB=NULL,beta=NULL,y.dbsum=NULL,ww.trisum=NULL,y.trisum=NULL,ww.trisum.d=NULL,y.trisum.d=NULL,bootstrap.iter=30){
  id <- as.character(id)
  ord <- order(match(id, unique(id)))
  id <- id[ord]
  Y <- as.matrix(Y)[ord, , drop = FALSE]
  X <- as.matrix(X)[ord, , drop = FALSE]
  Z <- as.matrix(Z)[ord, , drop = FALSE]
  if (is.null(q) == TRUE){
    q <- ncol(Z)
  }
  id <- factor(id, levels = unique(id))
  nis <- tabulate(id)
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z <- bdiag(lapply(Z_list, as.matrix))

  m <- length(nis)
  N <- sum(nis)
  if (is.null(d) == TRUE){
    d <- ncol(Y)
  }
  if (is.null(p) == TRUE){
    p <- ncol(X)
  }
  if (is.null(ww.trisum) == TRUE){
    ww.trisum <- triple_sum_ww(as.matrix(Z),nis,m,d=1,q)
  }
  if (is.null(ww.trisum.d) == TRUE){
    ww.trisum.d <- triple_sum_ww(as.matrix(Z),nis,m,d,q)
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


  if (ncol(Z) %% m != 0){
    stop("Error: The structure of Z incorrect!")
  }

  if (is.null(sigmaB)==TRUE){
    sigmaB <- Diagonal(d*q)
    #sigmaB <- matrix(0,d*q,d*q)
  }


  if (is.null(beta)==TRUE){
    beta <- solve(t(X) %*% X, t(X) %*% Y)
  }

  iter <- 0
  iterB <- NULL
  while (iter < max.iter) {
    iter <- iter + 1
    sigmaB.old <- sigmaB
    beta.old <- beta

    if (is.null(y.trisum.d)==TRUE){
      y.trisum.d <- triple_sum_y(as.matrix(Y),as.matrix(Z),X,as.matrix(beta),nis,m,d,q)
    }


    diag.vec <- diag(updateSigmaB(sigmaB,eta=NULL,ww.trisum.d,y.trisum.d,lambda.i=0,Y,Z,d,q,deltaB,max.iterB,tolB)$sigmaB)

    active.set.vec <- NULL
    active.set.mat <- NULL
    for (l in 1:d){
      Y_l <- Y[,l]
      beta_l <- beta[,l]
      sigmaB_l <- sigmaB[c((1+(l-1)*q):(l*q)),c((1+(l-1)*q):(l*q))]
      y.trisum <- triple_sum_y(as.matrix(Y_l),as.matrix(Z),X,as.matrix(beta_l),nis,m,d=1,q) #Could accelerate
      diag_l <- diag.vec[(1+(l-1)*q):(l*q)]
      weight_l <- 1/(diag_l^gamma.weight)
      lambda.weight <- lambdas[l] * diag(weight_l)
      sigmaB.sol <- updateSigmaB(sigmaB_l,eta,ww.trisum,y.trisum,lambda.weight,Y_l,Z,d=1,q,deltaB,max.iterB,tolB)
      sigmaB_l <- sigmaB.sol$sigmaB
      active.set.mat[[l]] <- which(diag(sigmaB_l) > threshold)
      active.set.vec <- c(active.set.vec,(active.set.mat[[l]]+q*(l-1)))
    }

    n <- length(active.set.vec)
    sigmaB <- estimateSigmaB(sigmaB=NULL,eta=NULL,Y,X,Z,beta,deltaB,active.set.vec,active.set.mat,n,nis,m,d,q)$sigmaB


    if (is.null(y.dbsum)==TRUE){
      y.dbsum <- double_sum_y(Y,X,beta,nis,m,d)
    }
    sigmaE <- updateSigmaE(y.dbsum,Z,sigmaB,nis,m,d,N,q,deltaE)

    beta <- updateBeta(Y,X,Z,tau,gamma,sigmaB,sigmaE,nis,d,p,q,N,bootstrap.iter)


    if (norm(sigmaB-sigmaB.old,"F")/max(1,norm(sigmaB.old,"F")) < tol){
      break
    }
  }


  return(list(SigmaB=sigmaB,SigmaE=sigmaE,beta=beta,iterB=iterB))
}
