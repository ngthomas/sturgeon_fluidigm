

#' compute probility of pairwise genotypes probabilities
#' 
#' From frequencies of different genotype categories and 
#' a genotyping error rate, compute probabilities of
#' unrelated and self-related pairs having certain pairs
#' of genotype categories.
#' @param p vector of frequencies of genotype categories
#' @param e scalar genotyping error rate
#' @return Returns a list with four components:
#' \itemize{
#'  \item \code{pu} a matrix of probabilities (rows are sample i's
#'  genotype category and columns for sample j's genotype category)
#'  given that i and j are unrelated.
#'  \item \code{ps} matrix of probs for case when i and j are samples from the same
#'  individual.
#'  \item \code{lr} log of the likelihood ratio of each genotype pair.  
#' }
gcat_pair_probs <- function(p, e) {
  ret <- list()
  ret$pu <- outer(p, p)
  ret$ps <- outer(p, p, function(g, h) ((1-e)^2) * g * (g == h) + (2*e - e^2) * g * h )
  lu <- log(ret$pu)
  ls <- log(ret$ps)
  ret$lr <- ls - lu
  ret
}


#' simulate pair genotype categories for a locus given input that is the output of gcat_pair_probs
#' 
#' Just simulates pairs that are self and unrelated and returns the logl ratio for each
#' @param a list that is like the output of gcat_pair_probs
#' @param reps the number of Monte Carlo samples to draw
#' @return Returns a list of two components.  U is a vector of length reps of logl-ratios
#' drawn according to their distirubtion given the pair is unrelated. S is a similar 
#' vector, but drawn according to their distribution given the two samples are from the
#' same individual.
sim_gcat_logls <- function(GC, reps) {
  list(
    U = sample(GC$lr, size = reps, replace = TRUE, prob = GC$pu),
    S = sample(GC$lr, size = reps, replace = TRUE, prob = GC$ps)
  )
}



#' do full sims from a list of genotype frequency category frequencies
#' 
#' @param GP a list of gc freqs
#' @param epsilon a genotyping error rate that is applied indisciminately to all loci
#' @param miss_rate  vector of rates of missing data by locus in same order as GP
#' @param reps
full_logl_sims <- function(GP, epsilon = 0.05, miss_rate, reps = 10^4, max_num_miss_loc = 15) { 
  GF <- lapply(GP, function(x) gcat_pair_probs(x, e = epsilon))
  tmp <- lapply(GF, function(x) sim_gcat_logls(x, reps = reps))
  
  Umat <- sapply(tmp, function(x) x$U)
  Smat <- sapply(tmp, function(x) x$S)
  
  Uvec <- rowSums(Umat)
  Svec <- rowSums(Smat)
  
  # now, get a matrix of reps rows, with, max_num_miss_loc columns, each row 
  # giving a sample of loci that would have been missing if just
  # one was missing, and if two were missing, and so forth
  missMat <- t(sapply(1:reps, function(y) sample(x = length(miss_rate), size = max_num_miss_loc, replace = FALSE, prob = miss_rate)))
  
  # now we are going to successively remove (NA them out) logls from the matrix of logls...
  ret <- list()
  ret[[1]] <- list(U = Uvec, S = Svec)
  for(i in 1:max_num_miss_loc) {
    Umat[ cbind(1:reps, missMat[,i])] <- NA
    Smat[ cbind(1:reps, missMat[,i])] <- NA
    
    ret[[i+1]] <- list(
      U = rowSums(Umat, na.rm = TRUE),
      S = rowSums(Smat, na.rm = TRUE)
    )
  }
  
  names(ret) <- seq(length(GP), length(GP) - max_num_miss_loc)
  ret
}