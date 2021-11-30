

## Functions for simulating data for testing


#' Simulate two related populations
#'
#' Simulates haplotypes from two populations that split 100 generations
#' ago, with 100 individuals from each population and 500 markers.
#'
#' @return A list containing two matrices of haplotypes.
#' @export
#'
#' @examples
simulate_populations <- function() {

  pop <- AlphaSimR::runMacs(nInd = 200,
                            nChr = 1,
                            split = 100,
                            segSites = 500)
  pop1 <- pop[1:100]
  pop2 <- pop[101:200]

  list(haplo1 = AlphaSimR::pullSegSiteHaplo(pop1),
       haplo2 = AlphaSimR::pullSegSiteHaplo(pop2),
       positions = pop@genMap[[1]][,1])
}
