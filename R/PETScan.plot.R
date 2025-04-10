#' Visualization of significant gene-peak pairs of score-based genome-wide association analysis
#' 
#' PETScan.plot() shows significant gene-peak pairs of score-based genome-wide 
#' association analysis using negative binomial models, 
#' with peak locations on the X-axis and gene locations on the Y-axis. 
#' Different levels of significance are indicated using different colors. 
#' The middle panel illustrates the number of genes associated with each peak, 
#' while the bottom panel displays the number of peaks linked to each gene.
#' 
#' @param geneID A vector of gene IDs from significant gene-peak pairs.
#' If geneID=k, its chromosome and position will be extracted from `eChr[k]` and `ePos[k]`, respectively.
#' @param peakID A vector of peak IDs from significant gene-peak pairs.
#' If peakID=k, its chromosome and position will be extracted from `pChr[k]` and `pPos[k]`, respectively.
#' @param scores A vector of p-values of each significant gene-peak pair.
#' @param scuts A vector indicating cutoffs of different significance levels.
#' @param cols A vector indicating colors for different significance levels.
#' @param eChr A vector specifying the chromosomes of genes, ordered by gene IDs.
#' @param ePos A vector specifying the genomic positions of the genes 
#' on their respective chromosomes, ordered by gene IDs.
#' @param pChr A vector specifying the chromosomes of the peaks, ordered by peak IDs.
#' @param pPos A vector specifying the genomic positions of the peaks 
#' on their respective chromosomes, ordered by peak IDs.
#' @param chroms A vector indicating the all the distinct chromosomes visualized.
#' @param plot.peak A logical indicating whether to count the number of peaks linked to each gene.
#' @param plot.gene A logical indicating whether to count the number of genes linked to each peak.
#' 
#' @references Sun, Wei. (2009). eQTL analysis by Linear Model. http://www.bios.unc.edu/~weisun/software/eMap.pdf
#' 
#' @examples
#' 
#' RNA = ENCODE_RNA
#' ATAC = ENCODE_ATAC
#' cvrt = ENCODE_cvrt[,2:4]
#' nb1 = PETScan(Ymat=RNA, X=cvrt, A=ATAC, perm=2, no_cores=8, seed=1) 
#' 
#' library(dplyr)
#' # Apply genomic control methods
#' nrow1 = nrow(ATAC)
#' score1 = nb1[1:nrow1, ]
#' gc1 = nb1[nrow1+1, ] / 0.456
#' gc1 = ifelse(gc1>1, gc1, 1)
#' score2 = t(t(score1) / gc1) %>% data.frame()
#' score2$ATAC = rownames(score2)
#' 
#' RNApos = ENCODE_RNApos
#' ATACpos = ENCODE_ATACpos
#' score3 = reshape2::melt(score2, id.vars="ATAC") %>%
#'   mutate(pvalue = 1 - pchisq(value, 1)) %>%
#'   mutate(padj = p.adjust(pvalue, method="BH")) %>%
#'   left_join(RNApos, by=c("variable"="gene")) %>%
#'   left_join(ATACpos, by=c("ATAC"="peak"))
#' 
#' PETScan.plot(score3$geneid, score3$peakid, score3$padj, 
#'           scuts = c(0.2, 0.1, 0.05, 0.01), 
#'           cols = c("yellow", "green", "blue", "red"),
#'           eChr = RNApos$g.chr, ePos = RNApos$g.start,
#'           pChr = ATACpos$p.chr, pPos = ATACpos$p.start, 
#'           chroms = c(1:22, "X"))
#' 
#' @export
PETScan.plot = function(geneID, peakID, scores, scuts, cols, 
                     eChr, ePos, pChr, pPos, chroms, plot.peak=TRUE, plot.gene=TRUE,
                     xlab="Peak Location", ylab="Gene Location"){
  
  if(length(geneID) == 0){
    stop("length(geneID)=0\n")
  }
  
  if(length(geneID) != length(peakID)){
    stop("length(geneID) != length(peakID)\n")
  }
  
  if(length(geneID) != length(scores)){
    stop("length(geneID) != length(scores)\n")
  }
  
  if(length(cols) < length(scuts)){
    stop("length(cols) < length(scuts)\n")
  }
  
  valid.chroms = c(1:90, "X", "Y")
  wrongChr = chroms[!(chroms %in% valid.chroms)]
  if(length(wrongChr)>0){
    stop(wrongChr, " are not valid chromosome labels\n")
  }
  
  todrop = which(!(eChr %in% chroms))
  eChr[todrop] = NA
  ePos[todrop] = NA
  eChr[eChr=="X"] = 99
  eChr[eChr=="Y"] = 100
  eChr = as.integer(eChr)
  ePos = as.numeric(ePos)
  
  todrop = which(!(pChr %in% chroms))
  pChr[todrop] = NA
  pPos[todrop] = NA
  pChr[pChr=="X"] = 99
  pChr[pChr=="Y"] = 100
  pChr = as.integer(pChr)
  pPos = as.numeric(pPos)
  
  chrs = union(unique(eChr), unique(pChr))
  chrs = sort(as.numeric(chrs))
  num.chrs = chrs[chrs <= 90]
  chr.chrs = c(num.chrs, "X", "Y")
  max.num = max(num.chrs)
  if(max.num > length(num.chrs)){
    str = "there is no peak/gene locaiton information for some chromosomes\n"
    stop(str)
  }
  eChr[eChr==99]  = max.num+1
  eChr[eChr==100] = max.num+2
  pChr[pChr==99]  = max.num+1
  pChr[pChr==100] = max.num+2
  
  chreMax = tapply(ePos, eChr, max, na.rm=TRUE, simplify = FALSE)
  chrmMax = tapply(pPos, pChr, max, na.rm=TRUE, simplify = FALSE)
  
  chrMax = numeric(length(chrs))
  for(i in 1:length(chrs)){
    ch = as.character(i)
    chrMax[i] = max(chreMax[[ch]], chrmMax[[ch]])
  }
  
  ek1 = list()
  scuts = sort(scuts, decreasing=TRUE)
  
  for(i in 1:(length(scuts)-1)){
    wkp = which(scores <= scuts[i] & scores > scuts[i+1])
    ek1[[i]] = data.frame(Gene_ID=geneID[wkp], Peak_ID=peakID[wkp])
  }
  i = length(scuts)
  wkp = which(scores <= scuts[i])
  ek1[[i]] = data.frame(Gene_ID=geneID[wkp], Peak_ID=peakID[wkp])
  
  if(plot.peak & plot.gene){
    layout(matrix(1:3,ncol=1), heights=c(8,2,2))
  }
  
  nChr = length(chrMax)
  chrLen = c(0, cumsum(chrMax))
  ep = ePos + chrLen[eChr]
  mp = pPos + chrLen[pChr]
  
  ymax = chrLen[nChr+1]
  bdr1 = -0.016*ymax
  bdr2 = -0.006*ymax
  
  par(mar=c(3,4,0,0))
  plot(c(bdr1,ymax*1.05), c(bdr1,ymax*1.05), type="n", xlab="",
       ylab="", main="", xaxt="n", yaxt="n", bty="n")
  mtext(xlab, side=1, line=1)
  mtext(ylab, side=2, line=1)
  for(i in 1:(length(scuts)-1)){
    gpos = ep[ek1[[i]]$Gene_ID]
    mpos = mp[ek1[[i]]$Peak_ID]
    points(mpos, gpos, col=cols[i], pch=21, cex=0.5)
  }
  
  i = length(scuts)
  gpos = ep[ek1[[i]]$Gene_ID]
  mpos = mp[ek1[[i]]$Peak_ID]
  points(mpos, gpos, col=cols[i], pch=21, cex=0.5)
  
  if(nChr>1){
    nchr.plot = floor(nChr/2)
    rect(rep(bdr1,nchr.plot), chrLen[seq(1,nChr,by=2)],
         rep(bdr2,nchr.plot), chrLen[seq(2,nChr+1,by=2)],
         border=NA, col="orange")
    rect(chrLen[seq(1,nChr,by=2)], rep(bdr1,nchr.plot),
         chrLen[seq(2,nChr+1,by=2)], rep(bdr2,nchr.plot),
         border=NA, col="orange")
  }
  
  kp = seq(1,nChr,by=2)
  ats = 0.5*(chrLen[-1] + chrLen[-length(chrLen)])
  mtext(chr.chrs[kp], at=ats[kp], side=1, line=-0.5, cex=0.8)
  mtext(chr.chrs[kp], at=ats[kp], side=2, line=-0.5, cex=0.8)
  
  lg = character(length(scuts)-1)
  for(i in 1:(length(scuts)-1)){
    lg[i] = sprintf("(%.1e, %.1e]", scuts[i+1], scuts[i])
  }
  i = length(scuts)
  lg[i] = sprintf("(0, %.1e]", scuts[i])
  
  legend(bdr1, ymax*1.1, lg, pch=21, col=cols,
         horiz = TRUE, bty="n", box.lty=2, cex=0.9)
  lines(c(0,ymax*1.01), rep(ymax*1.01,2), lty=2)
  lines(rep(ymax*1.01,2), c(0,ymax*1.01), lty=2)
  
  if(plot.peak){
    M.IDs = peakID[which(scores <= max(scuts))]
    et = table(M.IDs)
    mIDs = as.numeric(names(et))
    peak.cut = quantile(et, 0.95)
      
    par(mar=c(2,4,0,0))
    plot(mp[mIDs], et, xlim=c(bdr1,ymax*1.05),
         type="h", xlab="", ylab="Gene Count", xaxt="n", bty="n")
    abline(v=chrLen, lty=2, col="seagreen")
    mtext(chr.chrs[kp], at=ats[kp], side=1, line=0, cex=0.8)
    lines(c(0,ymax*1.01), rep(peak.cut,2), lty=2)
  }
  
  if(plot.gene){
    E.IDs = geneID[which(scores <= max(scuts))]
    mt = table(E.IDs)
    eIDs = as.numeric(names(mt))
    gene.cut = quantile(mt, 0.95)
    
    par(mar=c(2,4,0,0))
    plot(ep[eIDs], mt, xlim=c(bdr1,ymax*1.05),
         type="h", xlab="", ylab="Peak Count", xaxt="n", bty="n")
    abline(v=chrLen, lty=2, col="seagreen")
    mtext(chr.chrs[kp], at=ats[kp], side=1, line=0, cex=0.8)
    lines(c(0,ymax*1.01), rep(gene.cut,2), lty=2)
  }
}
