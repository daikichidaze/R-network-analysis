# 1.2 node list 
dg <- matrix(c(
  0,1,1,0,
  1,0,1,0,
  0,0,0,0,
  1,0,0,0),
  nrow = 4,
  ncol = 4,
  byrow = TRUE
  )
rownames(dg) <- letters[1:4]
colnames(dg) <- letters[1:4]
dg

# 1.3 edge list
e.list <- c(1,2,1,3,2,1,2,3,4,1)
e.matrix <- matrix(e.list, ncol = 2, byrow = TRUE)

# 1.4 file import
write.table(
  matrix(c(
    0,1,1,0,
    1,0,1,0,
    0,0,0,0,
    1,0,0,0),
    nrow = 4,
    ncol = 4,
    byrow = TRUE),
  file = "adj.txt",
  row.names = FALSE,
  col.names = FALSE
  )

g <- as.matrix(read.table("adj.txt"))

write.table(g,file = "adj.csv", row.names = FALSE, col.names = FALSE, sep=",")

g <- matrix(
  scan("adj.csv", sep=","),
  nrow = 4,
  ncol = 4,
  byrow = TRUE
)

# 1.7 bipartite graph
bg <- matrix(c(
  0,0,0,0,1,1,0,
  0,0,0,0,1,0,1,
  0,0,0,0,0,1,1,
  0,0,0,0,0,0,1,
  1,1,0,0,0,0,0,
  1,0,1,0,0,0,0,
  0,1,1,1,0,0,0),
  nrow=7)
rownames(bg) <- c("a", "b", "c", "d", "A", "B", "C")
colnames(bg) <- c("a", "b", "c", "d", "A", "B", "C")

bg %*% bg

# 1.8 sna
install.packages(("statnet"))

library(sna)
edgelist <- matrix(c(
  1,2,1,
  1,3,1,
  2,1,1,
  2,3,1,
  4,1,1),
  ncol = 3, byrow = TRUE
)

attr(edgelist, "n") <- 4
attr(edgelist, "vnames") <- letters[1:4]
edgelist

as.edgelist.sna(dg)
as.sociomatrix.sna(edgelist)

# multigraph
g1 <- matrix(c(
  0,1,0,
  1,0,1,
  0,1,0),
  nrow = 3)
g2 <- matrix(c(
  0,0,1,
  0,0,1,
  1,1,0),
  nrow = 3)
g <- array(dim = c(2,3,3))
g[1,,] <- g1
g[2,,] <- g2

g[2,,]

# weighted graph
wg <- matrix(c(
  0,2,0,4,
  2,0,3,1,
  0,3,0,0,
  4,1,0,0),
  nrow = 4,
  byrow = TRUE)

# 1.9 igraph
#install.packages("igraph")
library(igraph)

