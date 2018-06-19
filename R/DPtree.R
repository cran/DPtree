#' @importFrom Rdpack reprompt
NULL

#' Generating the standard D-P Tree prior
#'
#' \code{DPTreePrior} returns
#' a standard D-P Tree prior based on specified hyperparameters.
#'
#' @param m A positive integer. The finite approximation level for D-P tree. Default m=4.
#' @param z A positive number. On i-th level, the hyperparameter for D-P tree prior is
#' \eqn{z\times i^2}. Default z=1.
#' @return A list.
#' \item{a}{An array containing the hyperparameters of D-P trees.}
#' @examples
#' DPTreePrior(m=6, z=1)
#' @references
#' \insertRef{DPtree}{DPtree}
#' @export
DPTreePrior<-function(m=4, z=1){
  a <- array(NA, dim=c(2 ^ m, 2 ^ m, m))
  for(i in 1 : m){
    a[, ,i] <- z * i ^ 2
  }
  return(list(a=a))
}



#' Sampling a realized distribution from the D-P Tree.
#'
#' \code{RealizeDPTree} returns
#' a realized (copula) distribtuion sampled from the input D-P Tree.
#'
#' @param prior A list. Should be in same format as returned from  \code{DPTreePrior}.
#' @return An array of dimension  \eqn{2^m} by \eqn{2^m} by m. m is the approximation level.
#'  Realized Z's for all partitions at each level.
#' Three dimensions reprensent two marginals, and the level respectively.
#' @examples
#' RealizeDPTree(DPTreePrior(m=2, z=1))
#' @references
#' \insertRef{DPtree}{DPtree}
#' @export
RealizeDPTree <- function(prior){
  m <- dim(prior$a)[3]
  a <- prior$a
  Z <- array(NA, dim = c(2^m, 2^m, m))
  for (i in 1:m){                # the m+1-i level
    index <- ((0 : (2 ^ m - 1)) %/% 2 ^ (i - 1))
    for (j in (1 : 2 ^ (m - i) - 1) * 2){
      for (k in (1 : 2 ^ (m - i) - 1) * 2){
        alpha <- as.vector(a[2 ^ (i - 1) * (j + (1:2)), 2 ^ (i - 1) * (k + (1:2)), m + 1 - i])
        #print(alpha)
        d <- MCMCpack::rdirichlet(1, alpha)
        d <- matrix(d, 2, 2)
        Z[which(index == j), which(index == k), m + 1 - i] <- d[1, 1]
        #print(c(m+1-i,a0[1+2^(i-1)*j,1+2^(i-1)*k,m+1-i],a1[1+2^(i-1)*j,1+2^(i-1)*k,m+1-i]))
        Z[which(index == (j + 1)), which(index == (k + 1)), m + 1 - i]<-d[2, 2]
        Z[which(index == j), which(index == (k + 1)), m + 1 - i] <- d[1, 2]
        Z[which(index == j + 1), which(index == (k)), m + 1 - i] <- d[2, 1]
      }
    }
  }
  return(Z)
}

#' D-P tree posterior updating from a single copula observation.
#'
#' \code{DPTreePosterior} returns
#' the D-P tree posterior given input copula data.
#'
#' @param x An array of length 2. Single copula data observation.
#' Each element should be between 0 and 1.
#' @param prior A list. Should be in same format as returned from  \code{DPTreePrior}.
#' @param w A positive number. Weight of data for posterior updating. Default 1.
#' @return A list.
#' \item{a}{An array containing the hyperparameters of D-P trees.}
#' @examples
#' nsim = 1
#' rho = 0.9
#' data1 <- MASS::mvrnorm(n=nsim, mu=rep(0, 2), Sigma=matrix(c(1, rho, rho, 1), 2, 2))
#' data2 <- stats::pnorm(data1)
#' DPTreePosterior(x=data2, prior=DPTreePrior(m=4, z=1))
#' @references
#' \insertRef{DPtree}{DPtree}
#' @export
DPTreePosterior<-function(x, prior, w=1){
  m<-dim(prior$a)[3]
  a<-prior$a

  for(i in 1:m){
    index<-((0:(2^m-1))%/%2^(m-i))+1 #ith level
    #index.a<-((0:(2^m-1))%/%2^(m-i+1))+1
    cord<-c(min(trunc(x[1]*2^(i))+1,2^(i)),min(trunc(x[2]*2^(i))+1,2^(i)))

    a[which(index==cord[1]),which(index==cord[2]),i]<-a[which(index==cord[1]),which(index==cord[2]),i]+1*w
  }
  return(list(a=a))

}



#' D-P tree posterior updating from multiple copula observations.
#'
#' \code{DPTreePosteriorMulti} returns
#' the D-P tree posterior given input copula data.
#'
#' @param x An array of dimension n by 2. Multiple copula data observations,
#' with each row being a bivariate copula observation.
#' All elements should be between 0 and 1.
#' @param prior A list. Should be in same format as returned from  \code{DPTreePrior}.
#' @param w A positive number or an array of length n. Weight of data for posterior updating. Default 1.
#' @return A list.
#' \item{a}{An array containing the hyperparameters of D-P trees.}
#' @examples
#' nsim = 10
#' rho = 0.9
#' data1 <- MASS::mvrnorm(n=nsim, mu=rep(0, 2), Sigma=matrix(c(1, rho, rho, 1), 2, 2))
#' data2 <- stats::pnorm(data1)
#' DPTreePosteriorMulti(x=data2, prior=DPTreePrior(m=4, z=1))
#' @references
#' \insertRef{DPtree}{DPtree}
#' @export
DPTreePosteriorMulti<-function(x, prior, w=1){
  m <- dim(prior$a)[3]
  a <- prior$a

  for(i in 1:m){
    index <- ( (0:(2^m-1)) %/% 2 ^ (m - i)) + 1 #ith level
    #index.a<-((0:(2^m-1))%/%2^(m-i+1))+1
    ind1 <- (sapply(x[, 1], function(x){min(trunc(x * 2 ^ (i)) + 1, 2 ^ (i))}))
    ind2 <- (sapply(x[, 2], function(x){min(trunc(x * 2 ^ (i)) + 1, 2 ^ (i))}))
    cord <- cbind(ind1, ind2)
    cord1 <- do.call(rbind, apply(cord, 1,
                                  function(x){expand.grid(which(index == x[1]), which(index == x[2]))}))
    #print(cord1)
    cord2 <- as.matrix(plyr::count((cord1)))

    a[cbind(cord2[, 1:2], i)]<-a[cbind(cord2[, 1:2], i)] + cord2[, 3] * w
  }
  return(list(a=a))

}


#' Calculating sub-partition probabiltiy measures for a realized distribution from D-P tree.
#'
#' \code{DPTreeDensity} returns
#' the probablity measures in the finest sub-partitions of a realized distribution from D-P tree prior/posterior.
#'
#' @param Z An array of dimension of \eqn{2^m} by \eqn{2^m} by m, m being the approximation level.
#' Realized Z's for all partitions at each level,
#' as returned by \code{RealizeDPTree}.
#' @return A \eqn{2^m} by \eqn{2^m} matrix. Normalized measures for all \eqn{2^m} by \eqn{2^m} sub-partititons on copula space
#' given by the realized distribution from D-P tree.
#' @examples
#' dp.rlz <- RealizeDPTree(DPTreePrior(m=2, z=1))
#' DPTreeDensity(dp.rlz)
#' @references
#' \insertRef{DPtree}{DPtree}
#' @export
DPTreeDensity <- function(Z){
  m <- dim(Z)[3]
  d <- array(1, dim=c(dim(Z)[1], dim(Z)[2]))
  for (i in 1:m){
    d <- d * Z[, , i]
  }
  #print(sum(d))
  return(d / sum(d))
}


#' Calculating sub-partition probabiltiy measures for the posterior mean distribution from D-P tree.
#'
#' \code{DPTreePMeanDensity} returns
#' the probablity measures in the finest sub-partitions of the posterior mean from D-P tree.
#'
#' @param prior A list. D-P tree specification. Should be in same format as returned from \code{DPTreePrior} or \code{DPTreePosterior}.
#'
#' @return A \eqn{2^m} by \eqn{2^m} matrix. Normalized measures for all \eqn{2^m} by \eqn{2^m} sub-partititons on copula space
#' given by the posterior mean distribution from D-P tree.
#' @examples
#' DPTreePMeanDensity(DPTreePrior(m=2, z=1))
#' @references
#' \insertRef{DPtree}{DPtree}
#' @export
DPTreePMeanDensity<-function(prior){
  m <- dim(prior$a)[3]
  Z <- prior$a
  d <- array(1, dim=c(dim(Z)[1], dim(Z)[2]))
  n <- sum(prior$a[, , m] - m ^ 2)
  for (i in 1: (m - 1)){
    d <- d * Z[, , i] / (Z[, , i] - i ^ 2 + 4 * (i + 1) ^ 2)
  }
  d <- d * Z[, , m] / (n + 4 * 1 ^ 2)
  return(d)
}


#' The disitribution function for realized distribution from D-P tree.
#'
#' \code{pDPTreeRealize} returns
#' the value of distribution function of realized distribution from D-P tree at certain given point on copula space.
#'
#' @param d A \eqn{2^m} by \eqn{2^m} matrix, m being the approximating level. Normalized measures for all \eqn{2^m} by \eqn{2^m} sub-partititons on copula space
#' given by the realized distribution from D-P tree, as returned by \code{DPTreeDensity}.
#' @param x An array of dimension n by 2. The points on copula space for distribution function evluation. Should be between 0 and 1.
#' @return An array of length n. The values of CDF of the input D-P tree distribution evaluated at the input points.
#' @examples
#' pDPTreeRealize(DPTreePMeanDensity(DPTreePrior(m=2, z=1)),c(0.5,0.5))
#' @references
#' \insertRef{DPtree}{DPtree}
#' @export
pDPTreeRealize<-function(d, x){#x new data 2-dim; d density
  ff <- function(x){
    ngrid <- dim(d)[1]
    ind <- sapply(ceiling(x * ngrid), function(x){max(x, 1)})
    s <- numeric(4)
    #print(ind)
    s[1] <- ifelse(prod( (ind - 1)) == 0, 0, sum(d[1: (ind[1] - 1),1: (ind[2] - 1)]))
    s[2] <- ifelse( (ind[1] - 1) == 0, 0, sum(d[1: (ind[1] - 1), ind[2]]) * (x[2] - (ind[2] - 1) / ngrid)) * ngrid
    s[3] <- ifelse( (ind[2] - 1) == 0, 0, sum(d[ind[1], 1: (ind[2] - 1)]) * (x[1] - (ind[1] - 1) / ngrid)) * ngrid
    #print((x[1]-(ind[1]-1)/ngrid))
    s[4] <- (x[1] - (ind[1] - 1) / ngrid) * (x[2] - (ind[2] - 1) / ngrid) * d[ind[1], ind[2]] * ngrid ^ 2
    #print(s)
    return(sum(s))
  }
  F<-apply(matrix(x, ncol=2), 1, ff)
  return(F)
}


#' The disitribution function for realized distribution from D-P tree.
#'
#' \code{dDPTreeRealize} returns
#' the value of density function of realized distribution from D-P tree at certain given point on copula space.
#'
#' @param d A \eqn{2^m} by \eqn{2^m} matrix, m being the approximating level. Normalized measures for all \eqn{2^m} by \eqn{2^m} sub-partititons on copula space
#' given by the realized distribution from D-P tree, as returned by \code{DPTreeDensity}.
#' @param x An array of dimension n by 2. The points on copula space for density function evluation. Should be between 0 and 1.
#' @return An array of length n. The values of PDF of the input D-P tree distribution evaluated at the input points.
#' @examples
#' dDPTreeRealize(DPTreePMeanDensity(DPTreePrior(m=2, z=1)),c(0.5,0.5))
#' @references
#' \insertRef{DPtree}{DPtree}
#' @export
dDPTreeRealize <- function(d,x){#x new data 2-dim, can be matrix; d density
  ngrid<-dim(d)[1]
  ind<-(ceiling(x*ngrid))
  ind[ind==0]<-1
  ind<-matrix(ind,ncol=2)
  #print(ind)
  f<-apply(ind,1,function(y){d[y[1],y[2]]})*ngrid^2
  return(f)
}

#' Sample a copula observation from a realized distribution from D-P tree.
#'
#' \code{SampleDPTreeDensity} returns
#' a copula sample from a realized distribution from D-P tree.
#'
#' @param nsam A positive integer. The sample size.
#' @param d A \eqn{2^m} by \eqn{2^m} matrix, m being the approximating level.
#' Normalized measures for all \eqn{2^m} by \eqn{2^m} sub-partititons on copula space
#' given by the realized distribution from D-P tree, as returned by \code{DPTreeDensity}.
#'
#' @return An array of dimension nsam by 2. The values of PDF of the input D-P tree distribution evaluated at the input points.
#' @examples
#' SampleDPTreeDensity(10, DPTreePMeanDensity(DPTreePrior(m=2, z=1)))
#' @references
#' \insertRef{DPtree}{DPtree}
#' @export
SampleDPTreeDensity<-function(nsam, d){
  m <- log2(dim(d)[1])
  s <- sample(matrix(0: (2 ^ (2 * m) - 1), 2 ^ m, 2 ^ m), size=nsam, replace=T, prob=d)
  sample <- sapply(s, function(x){c(x %% (2 ^ m) + stats::runif(1), x %/% (2 ^ m) + stats::runif(1)) / 2 ^ m})
  return(t(sample))
}
