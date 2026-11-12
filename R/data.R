#' Simulated data to use in testing the colocmodern package
#'
#' @docType data
#'
#' @usage data(coloc_test_data)
#'
#' @format A list of four coloc-style datasets for a quantitative trait, D1 to
#'   D4, each over the same 500 SNPs, and `causals`, the position of each
#'   dataset's causal SNPs. D1 and D2 share a single causal variant. D3 and D4
#'   share one causal variant, and D3 has a second of its own. Use these as
#'   examples of what a coloc-style dataset for a quantitative trait should look
#'   like.
#'
#' @keywords datasets
#' @examples
#' data(coloc_test_data)
#' names(coloc_test_data)
#' str(coloc_test_data$D1)
#' check_dataset(coloc_test_data$D1) # should return NULL if data structure is ok
"coloc_test_data"
