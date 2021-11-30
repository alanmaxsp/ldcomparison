
#' Title
#'
#' @param filename_haps
#' @param filename_samples
#'
#' @return
#' @export
#'
#' @examples
read_shapeit <- function(filename_haps,
                         filename_samples) {


  ## transpose
}


#' Title
#'
#' @param haploA Vector of haplotypes for the first marker, coded 0/1.
#' @param haploB Vector of haplotypes for the second marker, coded 0/1.
#' @param stat What LD statistic to return ("D" for D or "r2" for r squared).
#'
#' @return Statistic (D or r squared).
#'
#' @examples
get_ld <- function(haploA,
                   haploB,
                   stat = "r2") {

  haploAB <- paste(haploA, haploB)
  haploAB[is.na(haploA)] <- NA
  haploAB[is.na(haploB)] <- NA

  nonmissingA <- sum(!is.na(haploA))
  nonmissingB <- sum(!is.na(haploB))
  nonmissingAB <- sum(!is.na(haploAB))

  pA <- sum(haploA, na.rm = TRUE)/nonmissingA
  pB <- sum(haploB, na.rm = TRUE)/nonmissingB

  pAB <- sum(haploAB == "1 1", na.rm = TRUE) / nonmissingAB

  D <- pAB - pA * pB

  if (stat == "r2") {

    r2 <- D^2 / (pA * (1 - pA) * pB * (1 - pB))

    result <- r2
  } else if (stat == "D") {
    result <- D
  }
  result
}



#' Title
#'
#' @param haplo Matrix of haplotypes, coded 0/1. Columns should be markers and
#' rows should be observed chromosomes.
#' @param stat What LD statistic to return ("D" for D or "r2" for r squared).
#'
#' @return A matrix of LD estimates.
#'
#' @examples
get_ld_matrix <- function(haplo,
                          stat = "r2") {

  fixed <- colSums(haplo) == 0 | colSums(haplo) == nrow(haplo)

  haplo <- haplo[,!fixed]

  n_qtl <- ncol(haplo)

  ld <- matrix(NA_real_,
               nrow = n_qtl,
               ncol = n_qtl)

  colnames(ld) <- colnames(haplo)
  rownames(ld) <- colnames(haplo)

  for (row_ix in 1:n_qtl) {
    for (col_ix in 1:row_ix) {
      ld[row_ix, col_ix] <- get_ld(haplo[, row_ix],
                                   haplo[, col_ix],
                                   stat = stat)
    }
  }

  ld
}



convert_ld_to_distance_pairs <- function(ld,
                                         positions,
                                         marker_names) {

  ld_df <- tibble::tibble(as.data.frame(ld))

  ld_df$marker1 <- colnames(ld_df)

  ld_long <- tidyr::pivot_longer(ld_df,
                                 -marker1,
                                 names_to = "marker2",
                                 values_to = "ld")

  ld_long <- na.exclude(ld_long)

  ld_long$marker1_position <- positions[match(ld_long$marker1, marker_names)]
  ld_long$marker2_position <- positions[match(ld_long$marker2, marker_names)]

  ld_long$distance <- abs(ld_long$marker1_position - ld_long$marker2_position)

  ld_long

}

get_phase_correlation <- function(pairwise_ld_population1,
                                  pairwise_ld_population2,
                                  positions) {
  cols <- c("marker1", "marker2", "ld")
  combined <- dplyr::inner_join(pairwise_ld_population1[, cols],
                                pairwise_ld_population2[, cols],
                                by = c("marker1", "marker2"))
  phase_cor <- cor(combined$ld.x, combined$ld.y)

  list(correlation = phase_cor,
       marker_pairs_used = nrow(combined))
}


group_ld_by_distance <- function(pairwise_ld,
                                 positions,
                                 window_size) {



}
