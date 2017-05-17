#' @title generate Jacobian matrix of different types 
#' @param N, number of nodes
#' @param C, connectance, 0 < C <= 1
#' @param adj_type, type of adjacency matrix
#' \describe{
#' \item{gnp}{ER random with proportion}
#' \item{gnm}{ER random with edge number}
#' \item{regular}{regular graph}
#' \item{fitness}{fitness model, scale-free graph}
#' }
#' @param mat_type, type of Jacobian matrix
#' \describe{
#' \item{pn}{positive-negative, antagonism}
#' \item{pp}{positive-positive, mutualism}
#' \item{nn}{negative-negative, competition}
#' \item{pz}{positive-zero, Commensalism}
#' \item{nz}{negative-zero, Amensalism}
#' \item{mixture}{mixture of the above five interaction-types}
#' \item{random}{random mixture}
#' }
#' @param dist_type, type of random distributions
#' \describe{
#' \item{normal}{normal distribution}
#' \item{uniform}{uniform distribution}
#' \item{binary}{binary distribution}
#' }
#' normal, uniform, bivalues
#' @param pn, proportion of positive-negative interactions
#' @param pp, proportion of positive-positive interactions
#' @param nn, proportion of negative-negative interactions
#' @param pz, proportion of positive-zero interactions
#' @param nz, proportion of negative-zero interactions
#' @return a Jacobian matrix
gen_Jacobian <- function(N, C, adj_type = c('gnp', 'gnm', 'regular', 'fitness'),
                         mat_type = c('random', 'pn', 'pp', 'nn', 'pz', 'nz', 'mixture'),
                         dist_type = c('normal', 'uniform', 'binary'), 
                         pn = 1, pp = 0, nn = 0, pz = 0, nz = 0,
                         mean = 0, sd = 1, ...) 
  {
  adj <- gen_adj(N, C, adj_type, ...)
  mat_type <- match.arg(mat_type)
  if (mat_type == 'pn') {  # positive-negative 
    pn = 1
    pp = nn = pz = nz = 0
  }
  else if (mat_type == 'pp') {  # positive-positive
    pp = 1
    pn = nn = pz = nz = 0
  }
  else if (mat_type == 'nn') {  # negative-negative
    nn = 1
    pn = pp = pz = nz = 0
  }
  else if (mat_type == 'pz') {  # positive-zero
    pz = 1
    pn = pp = nn = nz = 0
  }
  else if (mat_type == 'nz') {  # negative-zero
    nz = 1
    pn = pp = nn = pz = 0
  }
  else if (mat_type == 'random') {
    pn = 0.5
    pp = nn = 0.25
    pz = nz = 0
  }
  stopifnot(pn + pp + nn + pz + nz == 1)
  
  for (i in 2:N) {
    for (j in 1:(i-1)) {
      if (adj[i, j] == 1) { # if an interaction exists between i and j
        p = runif(1) # uniform random number in [0,1]
        if (p < pn / 2) {  # positive-negative - first half
          adj[i, j] = + gen_random(dist_type, mean, sd)
          adj[j, i] = - gen_random(dist_type, mean, sd)
        }
        else if (p < pn) {  # positive-negative - second half
          adj[i, j] = - gen_random(dist_type, mean, sd)
          adj[j, i] = + gen_random(dist_type, mean, sd)
        }
        else if (p < pn + pp) {  # positive-positive
          adj[i, j] = + gen_random(dist_type, mean, sd)
          adj[j, i] = + gen_random(dist_type, mean, sd)
        }
        else if (p < pn + pp + nn) {  # negative-negative
          adj[i, j] = - gen_random(dist_type, mean, sd)
          adj[j, i] = - gen_random(dist_type, mean, sd)
        }
        else if (p < pn + pp + nn + pz / 2) {  # positive-zero first half
          adj[i, j] = + gen_random(dist_type, mean, sd)
          adj[j, i] = 0
        }
        else if (p < pn + pp + nn + pz) {  # positive-zero second half
          adj[i, j] = 0
          adj[j, i] = + gen_random(dist_type, mean, sd)
        }
        else if (p < pn + pp + nn + pz + nz / 2) {  # negative-zero first half
          adj[i, j] = - gen_random(dist_type, mean, sd)
          adj[j, i] = 0
        }
        else if (p < pn + pp + nn + pz + nz) {  # negative-zero second half
          adj[i, j] = 0
          adj[j, i] = - gen_random(dist_type, mean, sd)
        }
        
      }
    }
  }
  adj
}

#' @title generate Adjacency matrix of different types 
gen_adj <- function(N, C, adj_type = c('gnp', 'gnm', 'regular', 'fitness'), ...) {
  require(igraph)
  adj_type <- match.arg(adj_type)
  if (adj_type == 'gnp') {
    G = sample_gnp(N, C)
  }
  else if (adj_type == 'gnm') {
    m = N * (N - 1) * C / 2  # number of (undirected) edges
    G = sample_gnm(N, m)
  }
  else if (adj_type == 'regular') {
    k = (N - 1) * C  # node (undirected) degree
    G = sample_k_regular(N, k)
  }
  else if (adj_type == 'fitness') {
    m = N * (N - 1) * C / 2  # number of (undirected) edges
    G = sample_fitness(m, (1:N)^-alpha) # alpha is the fitness of node
  }
  adj = as.matrix(as_adj(G))
}

#' @title generate a random number of different distribution types 
gen_random <- function(dist_type = c('normal', 'uniform', 'binary'), mean = 0, sd = 1) {
  dist_type <- match.arg(dist_type)
  if (dist_type == 'normal') {  # half normal distribution
    rand = abs(rnorm(1, mean = mean, sd = sd))
  }
  else if (dist_type == 'uniform') {  # uniform distribution
    rand = abs(runif(1, min = mean - sd, max = mean + sd))
  }
  else if (dist_type == 'binary') {
    rand = 1  # tricky!
  }
  rand
}