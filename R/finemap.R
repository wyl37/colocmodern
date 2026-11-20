## The coloc_finemap class: one trait's fine-mapping, as natural log Bayes
## factors for each signal and snp, which coloc() compares between two traits.

## The engines a coloc_finemap can come from. "external" is log Bayes factors
## given to as_coloc_finemap().
.finemap_engines <- c("wakefield", "susie", "external")

## Builds a coloc_finemap and checks it, the one place the class is given.
## A column named "null", SuSiE's column for no effect, is moved last, after the
## snps in their order. bf keeps no row names, as idx names its rows.
new_coloc_finemap <- function(bf, idx, snps, engine, args = list()) {
    if (is.matrix(bf)) {
        if (is.numeric(bf))
            storage.mode(bf) <- "double"
        rownames(bf) <- NULL
        nul <- colnames(bf) %in% "null"
        if (any(nul))
            bf <- bf[, c(which(!nul), which(nul)), drop = FALSE]
    }
    x <- validate_coloc_finemap(
        structure(
            list(bf = bf, idx = idx, snps = snps, engine = engine, args = args),
            class = "coloc_finemap"
        )
    )
    x$idx <- as.integer(x$idx)
    x
}

## Stops with a message for the first rule x breaks, and otherwise returns x.
## [[ reads each element by its exact name, where $ would take bfx for bf.
validate_coloc_finemap <- function(x) {
    fail <- function(...) stop(..., call. = FALSE)
    if (!is.list(x))
        fail("coloc_finemap must be a list; it is ", .what_is(unclass(x)))

    bf <- x[["bf"]]
    if (!is.matrix(bf) || !is.numeric(bf))
        fail("bf must be a numeric matrix of log Bayes factors, signals x snps; it is ",
             .what_is(bf))
    cn <- colnames(bf)
    if (ncol(bf) == sum(cn %in% "null"))
        fail("bf has no snps")
    if (is.null(cn))
        fail("bf must have snp ids as column names")
    if (any(is.na(cn) | cn == ""))
        fail("bf has an empty column name; name every column by snp id")
    if (any(!is.finite(bf)))
        fail("bf contains missing or infinite values")
    dup <- anyDuplicated(cn)
    if (dup > 0L)
        fail("snp id ", cn[dup], " is in bf twice")
    if ("null" %in% cn && cn[length(cn)] != "null")
        fail("the null column of bf must be the last")

    idx <- x[["idx"]]
    if (!is.numeric(idx) || !all(is.finite(idx)) || any(idx != round(idx)) ||
        any(idx < 1) || any(idx > .Machine$integer.max) || anyDuplicated(idx))
        fail("idx must be whole numbers above 0, with no missing value and no repeat")
    if (length(idx) != nrow(bf))
        fail("idx has ", length(idx), if (length(idx) == 1L) " value" else " values",
             " but bf has ", nrow(bf), if (nrow(bf) == 1L) " row" else " rows")

    snps <- x[["snps"]]
    if (!is.data.frame(snps))
        fail("snps must be a data.frame; it is ", .what_is(snps))
    if (!"snp" %in% names(snps))
        fail("snps must have a snp column")
    if (!is.character(snps[["snp"]]))
        fail("snps$snp must be a character vector of snp ids; it is ", .what_is(snps[["snp"]]))
    if (!identical(snps[["snp"]], cn[cn != "null"]))
        fail("snps$snp must be the column names of bf, in the same order, without null")

    engine <- x[["engine"]]
    if (!is.character(engine) || length(engine) != 1L || !engine %in% .finemap_engines)
        fail("engine must be one of: ", paste(.finemap_engines, collapse = ", "))
    args <- x[["args"]]
    if (!is.list(args) || is.object(args))
        fail("args must be a list; it is ", .what_is(args))
    x
}

##' Make a `coloc_finemap` from log Bayes factors
##'
##' A `coloc_finemap` holds one trait's fine-mapping in the form [coloc()]
##' compares: for each signal, the natural log Bayes factor of each snp against
##' no effect. [finemap()] makes one from a dataset; `as_coloc_finemap()` makes
##' one from natural log Bayes factors you have computed yourself.
##'
##' From a matrix, each row is a signal, `engine` is `"external"`, and `args`
##' is empty. A vector is treated as a single signal. A `coloc_finemap` is
##' returned as it is. If a `"null"` column (or element) is present, it is
##' placed last in `bf` and excluded from `snps`.
##'
##' @param x a numeric matrix of natural log Bayes factors (signals x snps),
##'   with snp ids as column names or given in `snps`; a numeric vector of
##'   natural log Bayes factors, named by snp id or with snp ids given in
##'   `snps`; or a `coloc_finemap`.
##' @param snps optional snp ids for the columns of `x`: a character vector, or
##'   a `data.frame` with a `snp` column. Other columns of a `data.frame` are
##'   kept in `snps`. When `x` has snp ids as names, `snps` must have the same
##'   ids in the same order, without `"null"`.
##' @param idx optional whole numbers above 0, one per signal, such as the
##'   SuSiE credible set each signal came from. By default, the row number.
##' @param ... not used. An argument given here stops with an error, so a
##'   misspelt `snps` or `idx` is not ignored.
##' @return A list of class `"coloc_finemap"` with these elements:
##'
##' * `bf`: a numeric matrix with one row per signal and one column per snp,
##'   named by snp id. A column named `"null"`, representing no effect (as given
##'   by SuSiE), comes last.
##' * `idx`: an integer vector with one entry per row of `bf`, such as the SuSiE
##'   credible set the row came from.
##' * `snps`: a `data.frame` with one row per snp, in the order of the columns
##'   of `bf`, without `"null"`. Its `snp` column holds the snp ids, and
##'   [finemap()] adds `position` and `z` when available.
##' * `engine`: `"wakefield"`, `"susie"`, or `"external"` (for log Bayes factors
##'   passed to `as_coloc_finemap()`).
##' * `args`: a list of the arguments the engine ran with.
##'
##' @seealso [finemap()], which makes a `coloc_finemap` from a dataset.
##' @export
##' @examples
##' bf <- matrix(c(0.1, 3.2, 0.4,
##'                2.5, 0.3, 6.1), nrow = 2, byrow = TRUE,
##'              dimnames = list(NULL, c("s1", "s2", "s3")))
##' as_coloc_finemap(bf)
##' as_coloc_finemap(c(s1 = 0.1, s2 = 3.2, s3 = 0.4))
##' ## the two rows came from SuSiE credible sets 2 and 5
##' as_coloc_finemap(unname(bf), snps = c("s1", "s2", "s3"), idx = c(2, 5))
as_coloc_finemap <- function(x, ...) {
    UseMethod("as_coloc_finemap")
}

##' @rdname as_coloc_finemap
##' @export
as_coloc_finemap.matrix <- function(x, snps = NULL, idx = seq_len(nrow(x)), ...) {
    .no_dots(...)
    if (is.null(snps)) {
        ## A matrix with no columns goes on to the validator, which says it has no snps
        if (ncol(x) > 0L && is.null(colnames(x)))
            stop("the matrix has no column names; name its columns by snp id, or give snps",
                 call. = FALSE)
        cn <- colnames(x)
        snps <- data.frame(snp = cn[!cn %in% "null"], stringsAsFactors = FALSE)
    } else {
        if (is.character(snps))
            snps <- data.frame(snp = snps, stringsAsFactors = FALSE)
        if (!is.data.frame(snps))
            stop("snps must be a character vector of snp ids or a data.frame with a snp column; it is ",
                 .what_is(snps), call. = FALSE)
        if (!"snp" %in% names(snps))
            stop("snps must have a snp column", call. = FALSE)
        ## R before 4.0 makes a factor of data.frame(snp = ids) by default
        if (is.factor(snps[["snp"]]))
            snps[["snp"]] <- as.character(snps[["snp"]])
        if (is.null(colnames(x))) {
            if (nrow(snps) != ncol(x))
                stop("snps has ", nrow(snps), " snps but x has ", ncol(x),
                     " log Bayes factors per signal", call. = FALSE)
            colnames(x) <- snps[["snp"]]
        }
    }
    new_coloc_finemap(x, idx, snps, "external")
}

##' @rdname as_coloc_finemap
##' @export
as_coloc_finemap.numeric <- function(x, snps = NULL, idx = 1L, ...) {
    .no_dots(...)
    if (is.null(names(x)) && is.null(snps))
        stop("the vector has no names; name it by snp id, or give snps", call. = FALSE)
    as_coloc_finemap(matrix(x, nrow = 1L, dimnames = list(NULL, names(x))),
                     snps = snps, idx = idx)
}

##' @export
as_coloc_finemap.coloc_finemap <- function(x, ...) {
    if (any(c("snps", "idx") %in% names(list(...))))
        stop("x is already a coloc_finemap; snps and idx apply only to a matrix or a vector",
             call. = FALSE)
    .no_dots(...)
    x
}

## Stops when ... holds anything, as a method's ... is there only because the
## generic has it, and a misspelt argument would otherwise be ignored
.no_dots <- function(...) {
    n <- ...length()
    if (n == 0L)
        return(invisible(NULL))
    nm <- names(list(...))
    if (is.null(nm))
        nm <- character(n)
    nm[nm == ""] <- "(unnamed)"
    stop("unused argument", if (n > 1L) "s", ": ", paste(nm, collapse = ", "), call. = FALSE)
}

##' @export
as_coloc_finemap.default <- function(x, ...) {
    stop("x is ", .what_is(x), "; as_coloc_finemap() takes a matrix of log Bayes factors, ",
         "a numeric vector or a coloc_finemap", call. = FALSE)
}

##' @export
print.coloc_finemap <- function(x, ...) {
    snps <- x[["snps"]][["snp"]]
    bf <- x[["bf"]]
    n <- nrow(bf)
    nsnps <- length(snps)
    cat("coloc finemap: ", x[["engine"]], ", ",
        n, if (n == 1L) " signal, " else " signals, ",
        nsnps, if (nsnps == 1L) " snp\n" else " snps\n", sep = "")
    if (n == 0L)
        cat("  no signals\n")
    for (i in seq_len(n)) {
        b <- bf[i, snps]
        top <- which.max(b)
        cat("  signal ", i, " (idx ", x[["idx"]][i], "): top snp ", snps[top],
            ", log BF ", format(signif(b[[top]], 3)), "\n", sep = "")
    }
    invisible(x)
}
