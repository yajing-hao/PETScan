#' Perform score-based genome-wide association analysis using negative binomial mixed models
#'
#' PETScan_Mixed() performs score-based genome-wide association analysis
#' across all RNA-Seq genes and ATAC-Seq peaks using negative binomial mixed models,
#' and returns a matrix where the (i,j) element of the matrix represents
#' the score test statistic for the i-th ATAC-Seq peak and the j-th RNA-Seq gene.
#'
#' @param Ymat A matrix of gene expression for RNA-Seq genes.
#' Rows and columns represent RNA-Seq gene names and samples, respectively.
#' @param X A matrix of covariates.
#' @param R A vector indicating the group of samples.
#' @param A A matrix of chromatin accessibility for ATAC-Seq peaks.
#' Rows and columns represent ATAC-Seq peak positions and samples, respectively.
#' @param perm An integer indicating the number of permutations for genomic control.
#' @param no_cores An integer indicating the number of cores for parallel calculations.
#' @param seed An integer specifying the seed for random number generation.
#' @param r0 A decimal specifying the threshold for reliability check.
#'
#' @return A matrix of score test statistics across all RNA-Seq genes and ATAC-Seq peaks.
#' Rows and columns represent ATAC-Seq peak positions and RNA-Seq gene names, respectively.
#' The last row contains the median of permuted score test statistics for each gene.
#'
#' @examples
#' RNA = GSE118165_RNA
#' ATAC = GSE118165_ATAC
#' cvrt = GSE118165_cvrt
#' nbm1 = PETScan_Mixed(Ymat=RNA, X=cvrt[,1:2], R=cvrt[,5], A=ATAC, perm=2, no_cores=8, seed=1)
#'
#' @import GLMMadaptive
#' @import parallel
#'
#' @export
PETScan_Mixed = function(Ymat, X=NULL, R, A, perm, no_cores=64, seed=1, r0=1) {
  library(parallel)
  R = factor(R)
  order0 = order(R)
  R = R[order0]
  Ymat = as.matrix(Ymat)[,order0]
  if (is.null(X)) {
    X = matrix(1, nrow=ncol(Ymat), ncol=1)
  } else {
    X = data.frame(X)[order0,]
    X = model.matrix(~., X)
  }
  A = as.matrix(A)[,order0]
  A = scale(t(A), center=TRUE, scale=FALSE)
  set.seed(seed)
  permA = function(A) {A[sample(1:nrow(A), replace=FALSE),]}
  Ap = do.call(cbind, replicate(perm, permA(A), simplify=FALSE))

  clust = makeCluster(no_cores)
  clusterEvalQ(clust, library(GLMMadaptive))
  res = parApply(clust, Ymat, 1, scorep.nbm, X, R, A, Ap, r0)
  stopCluster(clust)
  rownames(res) = c(unlist(colnames(A)), "GC")
  return(res)
}



#' Compute score test statistics for a given gene across peaks using negative binomial mixed models
#'
#' scorep.nbm() computes score test statistics for a given gene
#' across all ATAC-Seq peaks using negative binomial mixed models,
#' and return a vector where the i-th element of the vector represents
#' the score test statistic for the i-th ATAC-Seq peak and the given gene.
#'
#' @param Y A vector of gene expression for a given RNA-Seq gene.
#' @param X A matrix of covariates.
#' @param R A vector indicating the group of samples.
#' @param A A matrix of chromatin accessibility for ATAC-Seq peaks.
#' Rows and columns represent ATAC-Seq peak positions and samples, respectively.
#' @param Ap A matrix of permuted chromatin accessibility for ATAC-Seq peaks.
#' @param r0 A decimal specifying the filtering threshold.
#'
#' @return A vector of score test statistics for a given gene across all ATAC-Seq peaks
#' and the median of permuted score test statistics.
#'
#' @import GLMMadaptive
#'
#' @export
scorep.nbm = function(Y, X, R, A, Ap, r0) {
  error0 =  tryCatch({
    dat0 = data.frame(Y, R)
    nb1 = mixed_model(fixed = Y ~ 0 + X, random = ~ 1 | R, dat=dat0,
                      family=GLMMadaptive::negative.binomial, control=list(nAGQ=1))

    f0 = fitted(nb1, type="subject_specific")
    phi = exp(nb1[["phis"]])
    tau0 = nb1[["D"]][[1]]

    k0 = phi * f0 * (phi-f0) * (phi+Y) / (phi+f0)^3
    k1 = as.numeric( (phi * f0 * (phi+Y) / (phi+f0)^2) %*% model.matrix(~ factor(R) - 1) + 1/tau0 )
    k11 = rep(k1, times=table(R))
    k2 = as.numeric( ((2*phi * f0^2 + f0^2 * Y - phi * f0 * Y) / (phi + f0)^3) %*% model.matrix(~ factor(R) - 1) )
    k22 = rep(k2, times=table(R))

    ud1 = phi * (Y - f0) / (phi + f0)
    ud = ud1 - k0/k11/2
    U = colSums(ud * A)

    id1 = phi * f0 * (phi+Y) / (phi+f0)^2
    id2 = phi * f0 * (phi + Y) * (phi^2 - 4*phi*f0 + f0^2) / (phi + f0)^4
    id3 = k0 / k11
    id4 = id3 %*% t(id3)
    m0 = model.matrix(~ factor(R) - 1) %*% t(model.matrix(~ factor(R) - 1))
    id5 = id4*m0
    Dbb = diag(id1 + id2/k11/2) - id5/2
    ibb1 = t(X) %*% Dbb %*% A
    ibb2 = apply(A, 2, function(x) {t(x) %*% Dbb %*% x})

    Dbt = k0 / k11^2
    ibt = colSums(Dbt * A) / tau0

    d1 = f0 * (f0 - Y) / (phi + f0)^2
    d2 = (phi^2 * f0 * (4*f0-Y) + 2*phi* f0^2 * (2*Y-f0) - (f0^3*Y)) / (phi+f0)^4
    Dbp = d1 + d2/k11/2 - k0*k22/k11^2/2
    ibp = colSums(Dbp * A) * phi

    i1 = solve(nb1[["Hessian"]])
    i2 = rbind(ibb1, ibt, ibp)
    i3 = apply(i2, 2, function(x) {t(x) %*% i1 %*% x})
    I = ibb2 - i3
    r = i3 / ibb2
  }, error = function(e) {e})

  if ( inherits(error0, "error") ) {
    rep(NA, ncol(A)+1)
  } else if ( sum(is.na(r))>0 ) {
    rep(NA, ncol(A)+1)
  } else if ( any(r>r0) ) {
    rep(NA, ncol(A)+1)
  } else {
    Up = colSums(ud * Ap)
    ibb1p = t(X) %*% Dbb %*% Ap
    ibb2p = apply(Ap, 2, function(x) {t(x) %*% Dbb %*% x})
    ibtp = colSums(Dbt * Ap) / tau0
    ibpp = colSums(Dbp * Ap) * phi
    i2p = rbind(ibb1p, ibtp, ibpp)
    i3p = apply(i2p, 2, function(x) {t(x) %*% i1 %*% x})
    Ip = ibb2p - i3p
    rp = i3p / ibb2p
    Ip = ifelse(Ip>0 & rp<=r0, Ip, NA)

    t1 = U^2/I
    t2 = Up^2/Ip
    t2 = t2[!is.na(t2)]
    c(round(t1, 2), median(t2))
  }
}
