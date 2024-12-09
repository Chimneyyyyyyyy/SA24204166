#' use different model to cluster the graph
#'
#' @import fpc
#' @import kernlab
#' @param k number of clusters
#' @param cg1 complex matrix 1 generated from function DSBM
#' @param cg2 complex matrix 2 generated from function DSBM
#' @param g1 symmetric matrix generated from function DSBM
#' @param D diagonal matrix whose parameter is the degree of nodes
#' @param n nodes number
#'
#' @return list(vectors1,kmeansresult1,kmeansresult2)
#' @export
#'
#' @examples
#'  G<- DSBM(4,100,0.8,0.4)
#' different_model(4,G$cg1,G$cg2,G$g1,G$D,100)
#'
different_model<-function(k,cg1,cg2,g1,D,n) {
  A1 <- g1%*%t(g1)+t(g1)%*%g1
  A3<- solve(D)%*%A1
  D1 <- eigen(A3)
  P<-matrix(data=complex(real=rep(0,n*n),imaginary = rep(0,n*n)),nrow=n,ncol=n)
  A2<- eigen(solve(D)%*%cg1)
  for(i in 1:(k/2)){

    P<- P+as.matrix(A2$vectors[,i])%*%Conj(A2$vectors[,i])
  }
  km<- kmeans(cbind(Re(P),Im(P)),k)
  P<-matrix(data=complex(real=rep(0,n*n),imaginary = rep(0,n*n)),nrow=n,ncol=n)
  A2<- eigen(solve(D)%*%cg1)
  for(i in 1:(k/2)){
    P<- P+as.matrix((A2$vectors[,i])%*%Conj(t(A2$vectors[,i])))
  }
  km1<- kmeans(cbind(Re(P),Im(P)),k)
  P<-matrix(data=complex(real=rep(0,n*n),imaginary = rep(0,n*n)),nrow=n,ncol=n)
  A2<- eigen(solve(D)%*%cg2)
  for(i in 1:(k/2)){

    P<- P+as.matrix((A2$vectors[,i])%*%Conj(t(A2$vectors[,i])))
  }
  km2<- kmeans(cbind(Re(P),Im(P)),k)
  G<- list(kmeans(D1$vectors,k),km1,km2)
  names(G)<-c("DDSYMvector","km1","km2")
  return(G)

}


