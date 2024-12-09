#' Generate Directed SBM
#'
#' @import igraph
#' @importFrom stats kmeans
#' @importFrom stats runif
#' @param k number of cluster s
#' @param n nodes number
#' @param p probability of connecting nodes in the same cluster
#' @param q probability of connecting nodes in the different cluster
#'
#' @return list(symmetric matrix,complex matrix 1 ,complex matrix 2,diag(x))
#' @export
#'
#' @examples DSBM(4,100,0.8,0.4)
DSBM<-function(k,n,p,q){
  g <- make_empty_graph(n,directed=TRUE)
  g1<-matrix(c(0),nrow=n,ncol=n)
  realPart <- matrix(c(0),nrow=n,ncol=n)
  imagPart <- matrix(c(0),nrow=n,ncol=n)
  realPart1 <- matrix(c(0),nrow=n,ncol=n)
  imagPart1 <- matrix(c(0),nrow=n,ncol=n)
  realPart2 <- matrix(c(0),nrow=n,ncol=n)
  imagPart2 <- matrix(c(0),nrow=n,ncol=n)
  x1<-rep(0,n)
  # 平均分成十个类
  membership <- rep(1:n, each = n/10)
  # 同类中的节点以概论p连边
  intra_prob <- p

  # 不同类中的节点以概率q连边
  inter_prob <- q

  eta <- q
  # 添加边
  h <- 0
  for (i in 1:n) {
    for (j in i:n) {
      if (membership[i] == membership[j]&&i!=j) {
        if (runif(1) < intra_prob) {
          if(runif(1)<eta){
            g <- add_edges(g, c(i, j))
            g1[i,j]<- 1
            imagPart[i,j]<- 1
            imagPart[j,i]<- -1
            realPart1[i,j]<- cos(2*pi/n)
            imagPart1[i,j]<- sin(2*pi/n)
            realPart1[j,i]<- cos(2*pi/n)
            imagPart1[j,i]<- -sin(2*pi/n)
            x1[i]<- x1[i]+1
            x1[j]<- x1[j]+1
          }
          else{
            g <- add_edges(g, c(j, i))}
          g1[j,i]<- 1
          imagPart[i,j]<- -1
          imagPart[j,i]<- 1
          realPart1[i,j]<- cos(2*pi/n)
          imagPart1[i,j]<- -sin(2*pi/n)
          realPart1[j,i]<- cos(2*pi/n)
          imagPart1[j,i]<- sin(2*pi/n)
          x1[i]<- x1[i]+1
          x1[j]<- x1[j]+1
        }
      } else {
        if (runif(1) < inter_prob&&i!=j) {
          if(runif(1)<eta){
            g <- add_edges(g, c(i, j))
            g1[i,j]<- 1
            imagPart[i,j]<- 1
            imagPart[j,i]<- -1
            realPart1[i,j]<- cos(2*pi/n)
            imagPart1[i,j]<- sin(2*pi/n)
            realPart1[j,i]<- cos(2*pi/n)
            imagPart1[j,i]<- -sin(2*pi/n)
            x1[i]<- x1[i]+1
            x1[j]<- x1[j]+1
          }
          else{
            g <- add_edges(g, c(j, i))
            g1[j,i]<- 1
            imagPart[i,j]<- -1
            imagPart[j,i]<- 1
            realPart1[i,j]<- cos(2*pi/n)
            imagPart1[i,j]<- -sin(2*pi/n)
            realPart1[j,i]<- cos(2*pi/n)
            imagPart1[j,i]<- sin(2*pi/n)
            x1[i]<- x1[i]+1
            x1[j]<- x1[j]+1
          }
        }
      }
    }
  }
  complexMatrix1 <- matrix(data = complex(real=realPart,imaginary =imagPart),nrow=n,ncol=n)
  complexMatrix2 <- matrix(data = complex(real=realPart1,imaginary =imagPart1),nrow=n,ncol=n)
  D<- diag(x1+1)
  G<- list(g1,g,complexMatrix1,complexMatrix2,D)
  names(G)<-c("g1","g","cg1","cg2","D")
  return(G)
  }
