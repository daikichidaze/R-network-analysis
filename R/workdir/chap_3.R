# density
(2 * 6 ) / ( 5 * 4 )

6 / (5 *4)

Fig3.1 <- matrix(c(
  0,1,1,1,1,
  1,0,0,1,1,
  1,0,0,0,0,
  1,1,0,0,0,
  1,1,0,0,0),
  nrow = 5)

Fig3.2 <- matrix(c(
  0,1,1,1,1,
  0,0,0,1,1,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0),
  nrow=5,
  byrow=TRUE)

sum(Fig3.1)/(nrow(Fig3.1) * (nrow(Fig3.1) - 1))

sum(Fig3.2)/(nrow(Fig3.2) * (nrow(Fig3.2) - 1))

# transitivity
A <- matrix(c(
  0,1,1,1,1,1,0,0,
  1,0,1,0,0,0,1,0,
  1,1,0,0,0,0,0,0,
  1,0,0,0,1,0,0,1,
  1,0,0,1,0,1,0,0,
  1,0,0,0,1,0,0,0,
  0,1,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0),
  nrow=8)

A2 <- A %*% A
diag(A2) <- 0
A2
sum(A2 * A)/sum(A2)

# mutuality
Fig3.5 <- matrix(c(
  0,1,0,0,1,
  1,0,1,0,0,
  0,1,0,1,0,
  1,0,1,0,0,
  0,0,0,1,0),
  nrow=5,
  byrow=TRUE)

a <- Fig3.5 * t(Fig3.5)
acd <- Fig3.5 + t(Fig3.5)
acd[which(acd >= 1)] <- 1
sum(a)/sum(acd)

# sna
# density
library(sna)
gden(Fig3.1)
gden(Fig3.2)

# transitivity
gtrans(A)

# reciprocity
grecip(Fig3.5, measure = "dyadic.nonnull")

# mutuality
mutuality(Fig3.5)

grecip(Fig3.5, measure = "dyadic")

Fig3.8 <- matrix(c(
  0,0,1,1,0,
  0,0,0,1,1,
  0,0,0,1,0,
  0,0,0,0,0,
  0,1,0,0,0),
  nrow = 5, ncol = 5,
  byrow = TRUE)

sna::hierarchy(Fig3.8, measure = "krackhardt")
sna::efficiency(Fig3.8)
sna::lubness(Fig3.8)


#igraph
library(igraph)
g3.1 <- graph_from_adjacency_matrix(Fig3.1, mode = "undirected")
g3.2 <- graph_from_adjacency_matrix(Fig3.2)

igraph::edge_density(g3.1)
igraph::edge_density(g3.2)

g3.3 <- graph_from_adjacency_matrix(A, mode = "undirected")
igraph::transitivity(g3.3, type = "global")

g3.5 <- graph_from_adjacency_matrix(Fig3.5)
igraph::reciprocity(g3.5, mode = "ratio")

#3.7Example
library(statnet)
data(package = "sna")
data(package = "ergm")

?coleman

data(coleman)
coord1 <- gplot(coleman, g = 1)
#X11(width = 14, height = 7)
par(mfrow = c(1,2))
gplot(coleman, g= 1, coord = coord1, main = "Fall, 1957")
gplot(coleman, g=2, coord = coord1, main = "Spring, 1958")

sna::gden(coleman)
# 密度はそこそこ低いがわずかに上昇。
#「よく付き合う友達」という指標なので低いのかも
sna::gtrans(coleman) #推移性は下がっている。閉じた集団から開いた集団へ？
sna::grecip(coleman, measure="dyadic.nonnull") 
# 相互性は大きな変化はなし。微減は推移性と同様に、開いた集団へ？
sna::mutuality(coleman)
sna::connectedness(coleman)
sna::hierarchy(coleman, measure="krackhardt")
sna::efficiency(coleman)
sna::lubness(coleman)

# karate
?igraph::make_graph

(karate <- make_graph("Zachary"))
plot(karate)
edge_density(karate)
transitivity(karate)
