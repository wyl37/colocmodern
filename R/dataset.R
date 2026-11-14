## A coloc dataset constructor

## The elements with one value per snp, in the order check_dataset() compares
## their lengths and print() lists them
.per_snp_items <- c("beta", "varbeta", "pvalues", "MAF", "snp", "position", "z")

## The elements a coloc dataset may have. check_dataset() checks these and
## ignores anything else.
.dataset_items <- c(.per_snp_items, "N", "type", "s", "sdY", "LD")

##' Make a coloc dataset
##'
##' A coloc dataset holds one trait's summary statistics for one region, in
##' the form [finemap()] and [coloc()] read. `coloc_dataset()` makes one from
##' its elements, given as arguments, and `as_coloc_dataset()` makes one from a
##' list that holds them, such as a dataset made for the coloc package. Both
##' check the dataset with [check_dataset()], so a mistake is reported where
##' the dataset is made, not later inside [coloc()]. Make your datasets with
##' them before you call [finemap()], [coloc()] or [wakefield_bf()], which
##' stop when given a plain list.
##'
##' The result has class `"coloc_dataset"`, but is still a list with the same
##' elements, so code that reads `D$beta` works unchanged. `as_coloc_dataset()`
##' keeps elements coloc does not use.
##'
##' Every dataset needs `snp`. Which other elements it needs depends on the
##' engine that fine-maps it, and [check_dataset()] gives the details. In
##' brief, `"wakefield"` needs `type`, and `beta` and `varbeta` if you have
##' them. Without them, give `pvalues`, `MAF` and `N`, and `s` for a
##' case-control trait. A quantitative trait given as `beta` and `varbeta`
##' needs `sdY`, or `MAF` and `N` so coloc can estimate it. `"susie"` needs
##' `LD` and `N`, and either `z` or both `beta` and `varbeta`.
##'
##' A dataset is made with whichever elements you give, and their values are
##' checked. Whether it has what an engine needs is checked when the engine
##' runs, or by `check_dataset()` with that `engine`.
##'
##' If you change a dataset's elements after it is made, coloc's functions
##' check it again when they use it. `as_coloc_dataset()` returns a
##' coloc_dataset as it is, without checking it again; [check_dataset()]
##' checks one you have changed.
##'
##' A `type` given as a factor, or with a name, comes back as a plain string.
##'
##' @param snp character vector of snp ids
##' @param type `"quant"` for a quantitative trait, `"cc"` for
##'     case-control
##' @param beta,varbeta per-snp effect estimates, and their variances
##'     (std error squared)
##' @param pvalues per-snp p values
##' @param MAF per-snp minor allele frequencies
##' @param position per-snp positions
##' @param N number of samples
##' @param s for a case-control trait, the proportion of samples that
##'     are cases
##' @param sdY for a quantitative trait, the standard deviation of the
##'     trait
##' @param LD matrix of LD between snps, with snp ids as row and
##'     column names
##' @param z per-snp z-scores, optional: [finemap()] with `engine = "susie"`
##'     fine-maps from them in place of `beta / sqrt(varbeta)`
##' @param warn.minp warn if no p value is smaller than this
##' @return a list of class `"coloc_dataset"`
##' @seealso [check_dataset()] for what each element must look like
##' @export
##' @examples
##' D <- coloc_dataset(snp = paste0("s", 1:5),
##'                    type = "quant",
##'                    beta = c(0.1, 0.5, 0.2, 0, 0.05),
##'                    varbeta = rep(0.01, 5),
##'                    sdY = 1)
##' D
##'
##' ## a dataset made as a list for the coloc package
##' data(coloc_test_data)
##' D1 <- as_coloc_dataset(coloc_test_data$D1)
##'
##' ## mistakes are reported when the dataset is made
##' try(coloc_dataset(snp = paste0("s", 1:5),
##'                   type = "quant",
##'                   beta = c(0.1, 0.5, 0.2, NA, 0.05),
##'                   varbeta = rep(0.01, 5),
##'                   sdY = 1))
coloc_dataset <- function(snp = NULL, type = NULL, beta = NULL, varbeta = NULL,
                          pvalues = NULL, MAF = NULL,
                          position = NULL, N = NULL, s = NULL, sdY = NULL,
                          LD = NULL, z = NULL, warn.minp = 1e-6) {
    ## every argument defaults to NULL so that a missing snp or type gets
    ## check_dataset()'s message, not R's "argument is missing"
    d <- list(beta = beta, varbeta = varbeta, pvalues = pvalues, MAF = MAF,
              snp = snp, position = position,
              N = N, type = type, s = s, sdY = sdY, LD = LD, z = z)
    d <- d[!sapply(d, is.null)]
    new_coloc_dataset(d, warn.minp = warn.minp)
}

##' @rdname coloc_dataset
##' @param x a list with the elements of a coloc dataset, or a coloc_dataset
##' @export
as_coloc_dataset <- function(x, warn.minp = 1e-6) {
    if (inherits(x, "coloc_dataset")) return(x)
    new_coloc_dataset(x, warn.minp = warn.minp)
}

## Checks d and gives it the coloc_dataset class, the one place the class is
## given, so every coloc_dataset has been checked.
new_coloc_dataset <- function(d, warn.minp = 1e-6) {
    ## a data.frame and a susie fit are lists too, but neither is a dataset
    if (!is.list(d) || is.data.frame(d) || inherits(d, "susie"))
        .stop_not_a_dataset(d)
    d <- unclass(d)
    check_dataset(d, warn.minp = warn.minp)
    ## read by its label everywhere, so a factor or a name only gets in the way
    if (!is.null(d[["type"]]))
        d[["type"]] <- as.character(d[["type"]])
    structure(d, class = c("coloc_dataset", "list"))
}

##' @export
print.coloc_dataset <- function(x, ...) {
    ## the type by its label, as check_dataset() reads it: identical() would
    ## also compare a factor's or a named type's attributes
    cat("coloc dataset: ",
        if (!is.null(x[["type"]]))
            paste(if (isTRUE(x$type == "cc")) "case-control" else "quantitative", "trait, "),
        length(x$snp), " snps\n", sep = "")
    cat("  per-snp:", paste(intersect(setdiff(.per_snp_items, "snp"), names(x)), collapse = ", "), "\n")
    sc <- unlist(x[intersect(c("N", "s", "sdY"), names(x))])
    if (length(sc))
        cat("  ", paste(names(sc), signif(sc, 4), sep = " = ", collapse = ", "),
            "\n", sep = "")
    ## [[ matches the name exactly, where $ would take an element like LD_ref
    if (!is.null(x[["LD"]])) cat("  LD:", nrow(x[["LD"]]), "x", ncol(x[["LD"]]), "\n")
    other <- setdiff(names(x), .dataset_items)
    if (length(other)) cat("  other elements:", paste(other, collapse = ", "), "\n")
    invisible(x)
}
