A <- matrix(c(
  0,1,0,0,0,1,
  0,0,1,0,0,0,
  0,0,0,1,0,0,
  0,1,0,0,1,0,
  0,0,0,0,0,0,
  0,0,0,1,0,0),
  nrow = 6,
  byrow = TRUE
)

A %*% A

A %*% A %*% A

A %*% A %*% A %*% A

dijkstra <- function(matrix, vertex){
  L <- matrix
  L[which(L==0)] <- Inf
  diag(L) <- 0
  n <- nrow(matrix)
  d <- rep(Inf, n)
  d[vertex] <- 0
  M <- 1:n
  M <- M[-vertex]
  i <- vertex
  while (length(M) > 0) {
    for (j in 1:n){
      d[j] <- min(d[j], d[i] + L[i, j])
    }
    i <- M[which(d[M] == min(d[M]))[1]]
    M <- M[-which(M == i)]
  }
  return(d)
}

dijkstra(A, 1)

WF <- function(matrix) {
  L <- matrix
  L[which(L==0)] <- Inf
  diag(L) <- 0
  n <- nrow(matrix)
  for (k in 1:n){
    for (i in 1:n) {
      for (j in 1:n) {
        L[i,j] <- min(L[i,j], L[i,k] + L[k,j])
      }
    }
  }
  return(L)
}
WF(A)


install.packages(("statnet"))
library(sna)
geodist(A)

# weighted graph
wg <- matrix(c(
  0,2,0,4,
  2,0,3,1,
  0,3,0,0,
  4,1,0,0),
  nrow = 4,
  byrow = TRUE)
geodist(wg, ignore.eval = FALSE)

reachability(A)

is.connected(A, connected = "strong")

#igraph
install.packages(("devtools"))
require(devtools)

install_version('igraph', version='1.2.5', repos='https://cran.rstudio.org/')
require(igraph)

g1 <- graph_from_adjacency_matrix(A)
distances(g1, mode = "out")

g2 <- graph_from_adjacency_matrix(wg, weighted = TRUE)
distances(g2)

distances(g2, weights = c(2,4,2,3,1,3,4,1))
distances(g2, weights = NA)

shortest_paths(g1, from = 1, to = 5, mode = "out")$vpath

mean_distance(g1)          
mean_distance(g1, unconnected = FALSE)

distance_table(g1)

is_connected(g1, mode = "strong")

is_connected(g1, mode = "weak")
