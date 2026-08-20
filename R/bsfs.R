#' Convert one bSFS row into a blockwise pairwise-difference distribution
#'
#' Binomially splits the two samples' private-heterozygosity classes 50/50
#' and combines with the shared-heterozygosity and fixed-difference classes
#' to give, for one distinct (hetA, hetB, hetAB, fixed) combination observed
#' \code{count} times, the expected number of blocks contributing each
#' total pairwise-difference count. Direct port of \code{bsfsToS} from the
#' source Mathematica notebook.
#'
#' @param count number of blocks with this (hetA, hetB, hetAB, fixed) combination
#' @param hetA,hetB private heterozygous-site counts for each sample
#' @param hetAB shared heterozygous-site count
#' @param fixed fixed-difference count
#' @return numeric vector; element i is the expected block count for i - 1
#'   total pairwise differences contributed by this bSFS row
#' @export
bsfs_row_to_s <- function(count, hetA, hetB, hetAB, fixed) {
  het1 <- hetA + hetB
  tabhet1 <- dbinom(0:het1, size = het1, prob = 0.5)
  tabs1 <- tabhet1 * count
  target_len <- fixed + length(tabhet1)
  tabs1 <- c(rep(0, max(0, target_len - length(tabs1))), tabs1)  # pad on the left
  n <- length(tabs1) + hetAB
  right_pad <- c(tabs1, rep(0, n - length(tabs1)))
  left_pad  <- c(rep(0, n - length(tabs1)), tabs1)
  0.5 * right_pad + 0.5 * left_pad
}

#' Build a blockwise mutation-count distribution from a gIMble bSFS table
#'
#' Combines every distinct bSFS row for one species pair into a single
#' distribution over total pairwise differences per block - the direct
#' input format expected by \code{\link{fit_si}}, \code{\link{fit_im}},
#' \code{\link{fit_iim}}, \code{\link{fit_sc}} and
#' \code{\link{fit_demography}}. Direct port of \code{allbsfsToS}.
#'
#' @param bsfs_table a data.frame with columns \code{count}, \code{hetA},
#'   \code{hetB}, \code{hetAB}, \code{fixed} - one row per distinct bSFS
#'   tuple observed for a pair (i.e. the reshaped output of
#'   \code{gIMble query -b --bsfs})
#' @return a data.frame with columns \code{k} and \code{count}: the
#'   blockwise mutation-count distribution
#' @examples
#' toy <- data.frame(count = c(100, 50, 30), hetA = c(0, 1, 0),
#'                    hetB = c(0, 0, 0), hetAB = c(0, 0, 0), fixed = c(0, 0, 1))
#' bsfs_to_s_distribution(toy)
#' @export
bsfs_to_s_distribution <- function(bsfs_table) {
  rows <- Map(bsfs_row_to_s, bsfs_table$count, bsfs_table$hetA,
              bsfs_table$hetB, bsfs_table$hetAB, bsfs_table$fixed)
  max_len <- max(vapply(rows, length, integer(1)))
  padded <- lapply(rows, function(v) c(v, rep(0, max_len - length(v))))
  totals <- Reduce(`+`, padded)
  data.frame(k = seq_len(max_len) - 1L, count = totals)
}
