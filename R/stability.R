#' @title get stability measurements according to the Jacobian at equilibrium
#' @description stability measurements include:
#' 1. Rinf, asymptotic resilience, the negative of largest real eigenvalue
#' 2. R0, initial resilience(reactivity)
#' 3. Is, intrinsic stochastic invariability
#' 4. Id, intrinsic deterministic invariability. R0 < Is < Id < Rinf
#' @param J, the Jacobian at equilibrium
#' @return a list of stability measurements
#' @references Resilience, reactivity and variability : A mathematical comparison of ecological stability measures
get_stability <- function(J) {
  stopifnot(dim(J)[1] == dim(J)[2])
  n = dim(J)[1] # node number
  I = diag(1, n)
  Rinf = - max(Re(eigen(J)$values))
  R0 = - max(Re(eigen(J + t(J))$values)) / 2
  Is = 1 / (2 * norm(- solve(kronecker(I, J) + kronecker(J, I)), type = '2'))
  Id = 1 / norm(- solve(J), type = '2')
  list(Rinf = Rinf, R0 = R0, Is = Is, Id = Id)
}

#' @title localization of eigenvectors
get_localization <- function(eigenvectors) {
  stopifnot(dim(eigenvectors)[1] == dim(eigenvectors)[2])
  # assert that all eigenvectors is normalized
  if (! all(round(colSums(Mod(eigenvectors)^2),10) == 1))
    warning(round(colSums(Mod(eigenvectors)^2),10), 'Some eigenvector(s) is not normalized!')
  colSums(Mod(eigenvectors)^4)
}

#' @title get covariance matrix from Jacobian matrix
#' @param J Jacobian matrix
#' @param Sigma, the covariance matrix of environmental fluctuating
get_covariance <- function(J, Sigma = NULL) {
  n = dim(J)[1]
  if (is.null(Sigma)) {
    Sigma <- Matrix::Diagonal(n)
  }
  I = Matrix::Diagonal(n) #diag(1, n)
  - matrix(solve(kronecker(I, J) + kronecker(J, I)) %*% as.vector(Sigma), nrow = n, ncol = n)
}

#' @title get temporal stability from Covariance matrix
get_temporal_stability <- function(Cov) {
  Vc <- sum(Cov)
  Vs <- sum(sqrt(diag(Cov)))^2
  syn <- Vc / Vs
  list(Vc = Vc, Vs = Vs, syn = syn)
}

