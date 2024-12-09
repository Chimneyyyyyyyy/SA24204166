#' use Simpleherm model to cluster the graph
#'
#' @import fpc
#' @import kernlab
#' @param k number of clusters
#' @param cg2 complex matrix 2 generated from function DSBM
#' @param D diagonal matrix whose parameter is the degree of nodes
#' @param n nodes number
#'
#' @return list(vectors1,kmeansresult1,kmeansresult2)
#' @export
#'
#' @examples
#'  G<- DSBM(4,100,0.8,0.4)
#'  Simpleherm(4,G$cg2,G$D,100)
#'
Simpleherm<- function(k,cg2,D,n){
  P<-matrix(data=complex(real=rep(0,n*n),imaginary = rep(0,n*n)),nrow=n,ncol=n)
  A2<-diag(rep(1,n))-sqrt(solve(D))%*%cg2%*%sqrt(solve(D))
  D1<-eigen(A2)
  for(i in 1:(k/2)){

    P<- P+D1$vectors[,(n-k+1):n]%*%Conj(t(D1$vectors[,(n-k+1):n]))
  }
  km<- kmeans(cbind(Re(P),Im(P)),k)
  return(km)
}
#' use herm model to cluster the graph
#'
#' @import fpc
#' @import kernlab
#' @param k number of clusters
#' @param cg1 complex matrix 1 generated from function DSBM
#' @param D diagonal matrix whose parameter is the degree of nodes
#' @param n nodes number
#'
#' @return kmeans result
#' @export
#'
#' @examples
#'  G<- DSBM(4,100,0.8,0.4)
#'  HermRW(4,G$cg1,G$D,100)
#'
HermRW<- function(k,cg1,D,n){
  P<-matrix(data=complex(real=rep(0,n*n),imaginary = rep(0,n*n)),nrow=n,ncol=n)
  A2<-diag(rep(1,n))-sqrt(solve(D))%*%cg1%*%sqrt(solve(D))

  D1 <- eigen(A2,symmetric = TRUE)
  for(i in 1:(k/2)){

    P<- P+D1$vectors[,(n-k+1):n]%*%Conj(t(D1$vectors[,(n-k+1):n]))
  }

  km<- kmeans(cbind(Re(P),Im(P)),k)
  return(km)
}
#' use DDSYM model to cluster the directed graph
#'
#' @import fpc
#' @import kernlab
#' @param k number of clusters
#' @param g1 Adjacency Matrix  generated from function DSBM
#' @param n nodes number
#'
#' @return kmeans result
#' @export
#'
#' @examples
#'  G<- DSBM(4,100,0.8,0.4)
#'  DDSYM(4,G$g1,100)
#'
DDSYM<-function(k,g1,n) {
  x<-c(0)
  A1 <- g1%*%t(g1)+t(g1)%*%g1
  for(i in 1:n){
    x[i]<- sum(A1[i,])
  }
  D<-diag(x)
  A2<-diag(rep(1,100))-sqrt(solve(D))%*%A1%*%sqrt(solve(D))
  D1 <- eigen(A2,symmetric = TRUE)
  D2<-D1$vectors[,(n-k+1):n]
  return(kmeans(D2,k))

}
