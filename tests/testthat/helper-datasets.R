## Shared by the test files: a dataset in other forms.

## the same dataset as p-values, MAF and N
to_pvalues <- function(d) {
    d$pvalues <- 2 * pnorm(-abs(d$beta / sqrt(d$varbeta)))
    d[c("beta", "varbeta", "sdY")] <- NULL
    d
}

## the same dataset as case-control, and without position
as_cc <- function(d) {
    d$type <- "cc"
    d$s <- 0.3
    d$sdY <- NULL
    d
}
no_position <- function(d) {
    d$position <- NULL
    d
}
