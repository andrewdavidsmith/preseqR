# Copyright (C) 2016-2022 University of Southern California and
#                         Chao Deng and Andrew D. Smith and Timothy Daley
#
# Authors: Chao Deng
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

## predict the fraction of k-mers represented at least r times in the sample
kmer.frac <- function(n, r=2, mt=20) {
  return(preseqR.sample.cov(n=n, r=r-1, mt=mt))
}


#' @title Fraction of k-mers observed at least r times
#'
#' @description 'kmer.frac.curve' predicts the expected fraction of
#'     k-mers observed at least r times in a high-throughput
#'     sequencing experiment given the amount of sequencing
#'
#' @param n A two-column matrix. The first column is the frequency j
#'     = 1,2,...; and the second column is N_j, the number of k-mers
#'     observed exactly j times in the initial experiment. The first
#'     column must be sorted in an ascending order.
#' @param k The number of nucleotides in a k-mer.
#' @param read.len The average length of a read.
#' @param seq The amount of nucleotides sequenced.
#' @param r A positive integer. Default is 1.
#' @param mt An positive integer constraining possible rational
#'     function approximations. Default is 20.
#'
#' @return A two-column matrix. The first column is the amount of
#'     sequencing in an experiment. The second column is the estimate
#'     of the fraction of k-mers observed at least r times in the
#'     experiment.
#'
#' @details
#'
#' 'kmer.frac.curve' is mainly designed for metagenomics to evaluate
#' how saturated is a given metagenomic sequencing data set.
#'
#' 'kmer.frac.curve' is the fast version of
#' 'kmer.frac.curve.bootstrap'. The function does not provide the
#' confidence interval. To obtain the confidence interval along with
#' the estimates, one should use the function
#' 'kmer.frac.curve.bootstrap'.
#'
#' @author Chao Deng \email{chaodeng@usc.edu}
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1994). An introduction to the
#' bootstrap. CRC press.
#'
#' Deng C, Daley T, Calabrese P, Ren J & Smith AD (2020).
#' Predicting the Number of Bases to Attain Sufficient Coverage in
#' High-Throughput Sequencing Experiments.
#' Journal of Computational Biology, 27(7), 1130-1143
#'
#' @examples
#' # load library
#' library(preseqR)
#'
#' # import data
#' data(SRR061157_k31)
#'
#' # the fraction of 31-mers represented at least 10 times in an
#' # experiment when sequencing 1M, 10M, 100M, 1G, 10G, 100G, 1T
#' # nucleotides
#' kmer.frac.curve(n=SRR061157_k31, k=31, read.len=100,
#'                 seq=10^(6:12), r=10, mt=20)
#' @export
kmer.frac.curve <- function(n, k, read.len, seq, r=2, mt=20) {
  f <- kmer.frac(n, r=r, mt=mt)
  if (is.null(f))
    return(NULL)
  n[, 2] <- as.numeric(n[, 2])
  N <- n[, 1] %*% n[, 2]
  ## average number of k-mers per read
  m <- read.len - k + 1
  unit <- N / m * read.len
  # consistent vector-vector arithmetic
  unit <- as.numeric(unit)
  seq.effort <- seq / unit
  result <- matrix(c(seq, f(seq.effort)), ncol=2, byrow=FALSE)
  colnames(result) <- c("bases", paste("frac(X>=", r, ")", sep=""))
  return(result)
}


## predict the fraction of k-mers represented at least r times in the sample
kmer.frac.bootstrap <- function(n, r=2, mt=20, times=30, conf=0.95) {
  return(preseqR.sample.cov.bootstrap(n=n, r=r-1, mt=mt, times=times, conf=conf))
}

#' @title Fraction of k-mers observed at least r times with bootstrap
#'
#' @description Does the same thing as 'kmer.frac.curve’, predicting
#'     the expected fraction of k-mers observed at least r times in a
#'     high-throughput sequencing experiment given the amount of
#'     sequencing. But this function does it with bootstrapping.
#'
#' @param n A two-column matrix. The first column is the frequency j
#'     = 1,2,...; and the second column is N_j, the number of k-mers
#'     observed exactly j times in the initial experiment. The first
#'     column must be sorted in an ascending order.
#' @param k The number of nucleotides in a k-mer.
#' @param read.len The average length of a read.
#' @param seq The amount of nucleotides sequenced.
#' @param r A positive integer. Default is 1.
#' @param mt An positive integer constraining possible rational
#'     function approximations. Default is 20.
#' @param times The number of bootstrap samples.
#' @param conf The confidence level. Default is 0.95
#'
#' @return A four-column matrix. The first column is the amount of
#'     sequencing in an experiment. The second column is the estimate
#'     of the fraction of k-mers observed at least r times in the
#'     experiment. The third and fourth columns are the lower bounds
#'     and the upper bounds of the confidence intervals.
#'
#' @details This is the bootstrap version of ‘kmer.frac.curve’. The
#'     bootstrap sample is generated by randomly sampling the initial
#'     sample with replacement.  For each bootstrap sample, we
#'     construct an estimator. The median of estimates is used as the
#'     prediction for the number of species represented at least r
#'     times in a random sample.
#'
#'     The confidence interval is constructed based on a lognormal
#'     distribution.
#'
#' @author Chao Deng \email{chaodeng@usc.edu}
#'
#' @references
#'     Efron, B., & Tibshirani, R. J. (1994). An introduction to the
#'     bootstrap. CRC press.
#'
#'     Deng C, Daley T, Calabrese P, Ren J & Smith AD (2020).
#'     Predicting the Number of Bases to Attain Sufficient Coverage in
#'     High-Throughput Sequencing Experiments. Journal of
#'     Computational Biology, 27(7), 1130-1143
#'
#' @examples
#' # load library
#' library(preseqR)
#'
#' # import data
#' data(SRR061157_k31)
#'
#' # the fraction of 31-mers represented at least 10 times in an
#' # experiment when sequencing 1M, 10M, 100M, 1G, 10G, 100G, 1T
#' # nucleotides
#' kmer.frac.curve.bootstrap(n=SRR061157_k31, k=31, read.len=100,
#'                           seq=10^(6:12), r=10, mt=20)
#' @export
kmer.frac.curve.bootstrap <- function(n, k, read.len, seq, r=2, mt=20,
                                      times=30, conf=0.95)
{
  f <- kmer.frac.bootstrap(n, r=r, mt=mt, times=times, conf=conf)
  if (is.null(f))
    return(NULL)
  n[, 2] <- as.numeric(n[, 2])
  N <- n[, 1] %*% n[, 2]
  ## average number of k-mers per read
  m <- read.len - k + 1
  unit <- N / m * read.len
  # consistent vector-vector arithmetic
  unit <- as.numeric(unit)
  seq.effort <- seq / unit
  result <- matrix(c(seq, f$f(seq.effort), f$lb(seq.effort),
                     f$ub(seq.effort)), ncol=4, byrow=FALSE)
  colnames(result) <- c("bases", paste("frac(X>=", r, ")", sep=""),
                        "lb", "ub")
  return(result)
}
