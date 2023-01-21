D <- matrix(c(
  0,1,1,1,1,1,1,2,
  1,0,1,1,1,1,2,2,
  1,1,0,2,2,2,2,1,
  1,1,2,0,2,2,1,3,
  1,1,2,2,0,1,2,3,
  1,1,2,2,1,0,2,3,
  1,2,2,1,2,2,0,3,
  2,2,1,3,3,3,3,0),
  nrow = 8, ncol=8, byrow=TRUE)

1 / apply(D, 2, max)

1 / apply(D, 2, sum)

n <- nrow(D)
(n-1) / apply(D, 2, sum)

A <- matrix(c(
  0,1,1,1,1,1,1,0,
  1,0,1,1,1,1,0,0,
  1,1,0,0,0,0,0,1,
  1,1,0,0,0,0,1,0,
  1,1,0,0,0,1,0,0,
  1,1,0,0,1,0,0,0,
  1,0,0,1,0,0,0,0,
  0,0,1,0,0,0,0,0),
  nrow = 8, ncol=8, byrow=TRUE)
rowSums(A)

(evc <- abs(eigen(A)$vectors[,1]))

evc/max(evc)


# Page rank
B <- matrix(c(
  0,0,0,0,0,0,0,0,0,
  1,0,1,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,0,
  1,1,1,0,0,0,0,0,0,
  1,0,0,0,0,1,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,1,0,0,0),
  nrow=9,ncol=9,byrow=TRUE)

eigen(t(B))$values

B2 <- B
diag(B2)[which(rowSums(B2) == 0)] <- 1
(M <- t(B2 / rowSums(B2)))

t(B)

n <- nrow(M)
c <- 0.85
M <- (c * M) + ((1 - c) * matrix(1/n, nrow=n, ncol=n))
ev <- eigen(M)$vectors[,1]
rank <- ev / sum(ev)
as.numeric(rank)


1 / eigen(A)$values[1]
n <- nrow(A)
b <- 0.2
x <- rowSums(solve(diag(n) - b * A) %*% A)
a <- sqrt(n/sum(x^2))
a * x


# betweeenness
betweenness.centrality <- function(A){
  n <- nrow(A)
  Cb <- rep(0,n)
  for (s in 1:n){
    S <- c()
    P <- vector("list", n)
    g <- rep(0, n); g[s] <- 1
    d <- rep(-1, n); d[s] <- 0
    Q <- c()
    Q <- c(Q, s)
    while(length(Q) != 0){
      v <- Q[1]; Q <- Q[-1]
      S <- c(v, S)
      ws <- which(A[v,] == 1)
      if (length(Q) != 0){
        for (i in 1:length(ws)) {
          w <- ws[i]
          if (d[w] < 0){
            Q <- c(Q, w)
            d[w] <- d[v] + 1
          }
          if(d[w] == d[v] + 1){
            g[w] <- g[w] + g[v]
            P[[w]] <- c(P[[w]], v)
          }
        }
      }
    }
    b <- rep(0, n)
    while (length(S) != 0) {
      w <- S[1]; S <- S[-1]
      for (i in 1:length(P[[w]])) {
        v <- P[[w]][i]
        b[v] <- b[v] + (g[v] / g[w]) * (b[w] + 1)
      }
      if (w != s){
        Cb[w] <- Cb[w] + b[w]
      }
      
    }
  }
  
  Cb/2
}

Fig4.7 <- matrix(c(
  0,1,1,1,0,0,0,
  1,0,0,0,1,1,0,
  1,0,0,0,0,0,1,
  1,0,0,0,0,0,1,
  0,1,0,0,0,0,0,
  0,1,0,0,0,0,0,
  0,0,1,1,0,0,0),
  nrow = 7, ncol = 7, byrow = TRUE)

betweenness.centrality(Fig4.7)
