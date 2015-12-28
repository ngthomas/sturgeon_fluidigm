#' count up the number of pairs of discordant genotype category calls per locus
#' 
#' This takes as input a data frame that has, minimally, columns of 
#' \code{full.name}, \code{assay.name}, and \code{new.k}, though
#' these names can be set in the arguments to the function if they are named differently.
#' It counts up the number of pairs of non-missing genotype category calls that are 
#' discordant, as well as the total number of pairs, by assay.  
#' @param df a data frame holding the genotype category of each individual (one for each time
#' it was genotyped) at each assay.  Must have columns corresponding to the names in the 
#' remaining parameters.
#' @param nameCol quoted name of the column holding the individual name or ID
#' @param assayCol quoted name of the column holding the name of number of the assay
#' @param gcCol quoted name of the columns holding the genotype category of the individual. Note
#' that this will get coerced to a character, just to be safe about properly identifying missing
#' data.  
#' @return A list with two named components: \code{by_assay} is a data frame with columns named according to assayCol, 
#' and \code{discord_pairs} and \code{tot_pairs}. For each assay it gives the number of
#' discordant pairwise comparisons across chips and the total number of comparisons (over individuals).
#' The next component is \code{by_indiv} which gives the same, but summarized by individual over assays.
count_discordant_genotype_cats <- function(df,
                                           nameCol = "full.name",
                                           assayCol = "assay.name",
                                           gcCol = "new.k",
                                           missingCat = "0") {
  
  # this is a little function to count up the number of discordant pairs.  This
  # is sort of slow, but it works when there is more than just one pair (i.e. for
  # individuals that have been regenotyped multiple times)
  count_pairs <- function(x) {
    outer(x, x, "!=") %>%
      "["(lower.tri(.)) %>%
      sum
  }
  
  tmp <- df %>% 
    mutate_(., gcString = interp(~as.character(gcCol), gcCol = as.name(gcCol))) %>%  # make the geno cat a string.  (weird NSE here)
    filter(gcString != missingCat)  %>%   # filter out the missing data
    group_by_(assayCol, nameCol) %>%  # group to count up the number of genotypes of each individual at each assay
    mutate(numGenos = n())  %>%   # count it
    filter(numGenos > 1)  %>% # toss out those that don't have more than one replicate genotype
    summarise(num_discord = count_pairs(gcString), num_total = n() * (n() - 1) / 2) %>%
    ungroup
  
  list(by_assay = tmp %>% 
         group_by_(assayCol) %>%
         summarise(discord_pairs = sum(num_discord), tot_pairs = sum(num_total)),
       by_indiv = tmp %>% 
         group_by_(nameCol) %>%
         summarise(discord_pairs = sum(num_discord), tot_pairs = sum(num_total))
       )
    
}