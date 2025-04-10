#' Perform score-based genome-wide association analysis using negative binomial models
#'
#' PETScan() performs score-based genome-wide association analysis
#' across all RNA-Seq genes and ATAC-Seq peaks using negative binomial models,
#' and returns a matrix where the (i,j) element of the matrix represents
#' the score test statistic for the i-th ATAC-Seq peak and the j-th RNA-Seq gene.
#'
#' @param Ymat A matrix of gene expression for RNA-Seq genes.
#' Rows and columns represent RNA-Seq gene names and samples, respectively.
#' @param X A matrix of covariates.
#' @param A A matrix of chromatin accessibility for ATAC-Seq peaks.
#' Rows and columns represent ATAC-Seq peak positions and samples, respectively.
#' @param perm An integer indicating the number of permutations for genomic control.
#' @param no_cores An integer indicating the number of cores for parallel calculations.
#' @param seed An integer specifying the seed for random number generation.
#'
#' @return A matrix of score test statistics across all RNA-Seq genes and ATAC-Seq peaks.
#' Rows and columns represent ATAC-Seq peak positions and RNA-Seq gene names, respectively.
#' The last row contains the median of permuted score test statistics for each gene.
#'
#' @examples
#' RNA = ENCODE_RNA
#' ATAC = ENCODE_ATAC
#' cvrt = ENCODE_cvrt[,2:4]
#' nb1 = PETScan(Ymat=RNA, X=cvrt, A=ATAC, perm=2, no_cores=8, seed=1)
#'
#' @import MASS
#' @import parallel
#'
#' @export
PETScan = function(Ymat, X=NULL, A, perm, no_cores=64, seed=1) {
  library(parallel)
  Ymat = as.matrix(Ymat)
  if (is.null(X)) {
    X = matrix(1, nrow=ncol(Ymat), ncol=1)
  } else {
    X = model.matrix(~., data.frame(X))
  }
  A = as.matrix(A)
  set.seed(seed)
  permA = function(A) {A[,sample(1:ncol(A), replace=FALSE)]}
  Ap = do.call(rbind, replicate(perm, permA(A), simplify=FALSE))

  clust = makeCluster(no_cores)
  clusterEvalQ(clust, library(MASS))
  res = parApply(clust, Ymat, 1, scorep.nb, X, A, Ap)
  stopCluster(clust)
  rownames(res) = c(unlist(rownames(A)), "GC")
  return(res)
}



#' Compute score test statistics for a given gene across peaks using negative binomial models
#'
#' scorep.nb() computes score test statistics for a given gene
#' across all ATAC-Seq peaks using negative binomial models,
#' and returns a vector where the i-th element of the vector represents
#' the score test statistic for the i-th ATAC-Seq peak and the given gene.
#'
#' @param Y A vector of gene expression for a given RNA-Seq gene.
#' @param X A matrix of covariates.
#' @param A A matrix of chromatin accessibility for ATAC-Seq peaks.
#' Rows and columns represent ATAC peak positions and samples, respectively.
#' @param Ap A matrix of permuted chromatin accessibility for ATAC-Seq peaks.
#'
#' @return A vector of score test statistics for a given gene across all ATAC-Seq peaks
#' and the median of permuted score test statistics.
#'
#' @import MASS
#'
#' @export
scorep.nb = function(Y, X, A, Ap) {
  error0 =  tryCatch({
    nb1 = glm.nb(Y ~ 0 + X)
  }, error = function(e) {e})

  if (inherits(error0, "error")) {
    return(rep(NA, nrow(A)+1))
  } else {
    fitted0 = nb1$fitted.values
    theta0 = nb1$theta

    ud = theta0 * (Y - fitted0) / (theta0 + fitted0)
    id = theta0 * fitted0 * (theta0 + Y) / (theta0 + fitted0)^2
    Xd = id * X
    IM = diag(length(Y)) - X %*% solve(crossprod(X, Xd), t(Xd))
    E2 = IM %*% t(A)
    U = colSums(ud * E2)
    I = colSums(id * E2^2)

    E2p = IM %*% t(Ap)
    Up = colSums(ud * E2p)
    Ip = colSums(id * E2p^2)
    tp = Up^2/Ip
    tp = tp[!is.na(tp)]

    return(c(round(U^2/I, 2), median(tp)))
  }
}
