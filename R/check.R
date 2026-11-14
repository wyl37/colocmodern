##' Check a coloc dataset
##'
##' `check_dataset()` checks that coloc can use a dataset. It stops at the
##' first problem it finds, with a message that names the element and says
##' what is wrong. You seldom need to call it, as [coloc_dataset()] and
##' [as_coloc_dataset()] check a dataset when they make it, and [coloc()],
##' [finemap()] and [wakefield_bf()] check every dataset they are given. Call
##' it yourself to check your datasets before you analyse them, for example
##' every region of a study before a loop over the regions, or to check a
##' dataset you have changed.
##'
##' A coloc dataset is a list that describes one trait in one region. Some of
##' its elements have one value per snp, and these must all have the same
##' length:
##'
##' * `snp`, the snp ids, which coloc matches the two traits on;
##' * `beta`, the regression coefficient of each snp, and `varbeta`, its
##'   variance, the standard error squared;
##' * `pvalues`, the p-value of each snp, and `MAF`, its minor allele
##'   frequency;
##' * `position`, the position of each snp;
##' * `z`, the z-score of each snp, which [finemap()] with `engine = "susie"`
##'   fine-maps from in place of `beta / sqrt(varbeta)`.
##'
##' The other elements are single values that describe the study, and the LD
##' matrix:
##'
##' * `type`, `"quant"` for a quantitative trait or `"cc"` for a case-control
##'   trait;
##' * `N`, the number of samples;
##' * `s`, for a case-control trait, the proportion of samples that are cases;
##' * `sdY`, for a quantitative trait, the standard deviation of the trait;
##' * `LD`, the correlation matrix of the snps, with snp ids as its row and
##'   column names.
##'
##' Every dataset needs `snp`. Which other elements it needs depends on the
##' engine that fine-maps it:
##'
##' * `"wakefield"` needs `type`. Give `beta` and `varbeta` if you have them.
##'   Without them, it approximates them from `pvalues`, `MAF` and `N`, and a
##'   case-control trait also needs `s`. A quantitative trait given as `beta`
##'   and `varbeta` needs `sdY`, or `MAF` and `N` so coloc can estimate `sdY`.
##' * `"susie"` needs `LD` and `N`, and either `z` or both `beta` and
##'   `varbeta`.
##'
##' `check_dataset()` first checks the dataset's structure: it is a list and
##' not a data.frame, `snp` is given, with each snp id once and not as a
##' factor, the elements `req` names are given, and the per-snp elements have
##' the same length.
##'
##' It then checks the values of some elements. With `engine = "none"`,
##' these are all the elements that are given, apart from `z`. With
##' `"wakefield"` or `"susie"`, they are the elements that engine reads, which
##' must all be given. So a check for `"wakefield"` leaves `LD` alone, and
##' checking `LD` takes far longer than checking the rest of a dataset. An
##' element whose values are checked must have no missing values, and:
##'
##' * `type` must be one value, `"quant"` or `"cc"`;
##' * `s` must be one number between 0 and 1;
##' * `MAF` must be numbers above 0 and below 1;
##' * `N` must be one finite number above 1;
##' * `beta` and `varbeta` must be numeric, with no infinite values, and
##'   `varbeta` must be above 0;
##' * `pvalues` must be numbers above 0 and at most 1;
##' * `sdY` must be one finite number above 0;
##' * `position` must hold finite whole numbers above 0;
##' * `LD` must be a numeric matrix and square, with row names that match its
##'   column names, in the same order, and cover every snp. For the dataset's
##'   snps it must be a correlation matrix, with values from -1 to 1,
##'   symmetric, and 1 on the diagonal;
##' * `z`, which only `"susie"` reads, must be a finite number for each snp.
##'
##' It also warns, without stopping, when no p-value is below `warn.minp`, as
##' the trait may then have no association in the region. A [coloc_dataset()]
##' gets this warning once, when it is made.
##'
##' @param d a dataset, as a list or a [coloc_dataset()], but not a data.frame
##' @param suffix a label for the dataset, such as 1 or 2, which begins each
##'     message as "dataset 1: ". With "" messages have no label
##' @param engine the engine the dataset is checked for, `"wakefield"` or
##'     `"susie"`, or `"none"`, the default, for no engine
##' @param req other elements that must be present, such as `"position"`
##' @param warn.minp warn when no p-value is below this
##' @return `NULL` when it finds no problem. Otherwise it stops with an error
##' @importFrom stats pnorm
##' @export
##' @author Chris Wallace
##' @examples
##' data(coloc_test_data)
##' check_dataset(coloc_test_data$D1)     # NULL: no problems found
##' check_dataset(coloc_test_data$D1, engine = "susie")
##' bad <- coloc_test_data$D1
##' bad$beta[1] <- NA
##' try(check_dataset(bad, 1))
check_dataset <- function(d,
                          suffix="",
                          engine=names(.engine_reads),
                          req=NULL,
                          warn.minp=1e-6) {
    ## A label, such as 1 or 2, starts each message as "dataset 1: "; "", NULL,
    ## NA or more than one value is no label. No call is shown, as the message
    ## names the dataset.
    dataset <- if (length(suffix) == 1L && !is.na(suffix) && !identical(as.character(suffix), ""))
                   suffix
    prefix <- if (is.null(dataset)) "" else paste0("dataset ", dataset, ": ")
    fail <- function(...) stop(prefix, ..., call. = FALSE)

    engine <- tryCatch(match.arg(engine), error = function(e)
        stop("engine must be one of: ", paste(names(.engine_reads), collapse = ", "), call. = FALSE))
    if (!is.list(d) || is.data.frame(d)) # a data.frame repeats type, sdY and N on every row
        .stop_not_a_dataset(d, dataset)
    nd <- intersect(names(d), .dataset_items)          # the elements coloc reads

    ## The checks stop at the first problem, in two blocks: the dataset's
    ## structure, then the values of the elements in `read`.

    ## 1. The structure, on every call
    absent <- setdiff(c("snp", req), nd)
    if (length(absent))
        fail("missing required element(s) ", paste(absent, collapse = ", "))
    if (any(duplicated(d$snp)))
        fail("duplicated snps found")
    if (is.factor(d$snp))
        fail("snp should be a character vector but is a factor")
    ## each per-snp element but z has as many values as snp
    for (v in intersect(setdiff(.per_snp_items, "z"), nd))
        if (length(d[[v]]) != length(d$snp))
            fail("snp has ", length(d$snp), if (length(d$snp) == 1L) " value" else " values",
                 " but ", v, " has ", length(d[[v]]))

    ## 2. The values of the elements in `read`, which are what the engine
    ## reads. The engine's entry stops when an element it needs is absent.
    read <- c("snp", .engine_reads[[engine]](d, nd, fail))
    reads <- function(v) v %in% read
    ## no missing values - make people clean their own data rather than make
    ## assumptions here for datasets I don't know. z has its own rule, below
    for (v in setdiff(intersect(nd, read), "z"))
        if (anyNA(d[[v]]))
            fail(v, " contains missing values")
    if (reads("type") && !(length(d$type) == 1L && d$type %in% c("quant", "cc")))
        fail("type must be quant or cc")
    if (reads("s")) {
        if (length(d$s) != 1L)
            fail("s must be a single number between 0 and 1")
        if (!is.numeric(d$s) || d$s <= 0 || d$s >= 1)
            fail("s must be between 0 and 1")
    }
    if (reads("MAF") && (!is.numeric(d$MAF) || any(d$MAF <= 0 | d$MAF >= 1)))
        fail("MAF should be a numeric, strictly >0 & <1")
    if (reads("N")) {
        if (!is.numeric(d$N) || length(d$N) != 1L)
            fail("N must be a single number; coloc does not take per-snp sample sizes")
        if (!is.finite(d$N) || d$N <= 1)
            fail("N must be finite and > 1")
    }
    for (v in intersect(c("beta", "varbeta"), read))
        if (any(is.infinite(d[[v]])))
            fail(v, " has infinite values")
    ## a varbeta that isn't a number is refused below
    if (reads("varbeta") && is.numeric(d$varbeta) && any(d$varbeta == 0))
        fail("varbeta must be > 0")
    if (reads("position") &&
        !(is.numeric(d$position) && all(is.finite(d$position)) && all(d$position > 0) &&
          all(d$position == round(d$position))))
        fail("position must hold finite whole numbers above 0")
    if (reads("LD")) {
        if (length(dim(d$LD)) < 2L)
            fail("LD must be a matrix")
        if (nrow(d$LD) != ncol(d$LD))
            fail("LD must be square, but it has ", nrow(d$LD), ngettext(nrow(d$LD), " row", " rows"),
                 " and ", ncol(d$LD), ngettext(ncol(d$LD), " column", " columns"))
        if (!identical(colnames(d$LD), rownames(d$LD)))
            fail("LD rownames != colnames")
        if (length(setdiff(d$snp, colnames(d$LD))))
            fail("LD must cover every snp")
        if (length(dim(d$LD)) > 2L)
            fail("LD must be a matrix")
        if (!(is.matrix(d$LD) && is.numeric(d$LD)))
            fail("LD must be a numeric matrix")
        ## only the dataset's snps, copied out only when LD covers more; with
        ## as many, they are the same snps
        LD <- if (ncol(d$LD) > length(d$snp)) d$LD[d$snp, d$snp] else d$LD
        if (!.is_correlation(LD))
            fail("LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal")
    }
    if (reads("pvalues")) {
        if (!is.numeric(d$pvalues))
            fail("pvalues must be numeric")
        if (any(d$pvalues <= 0 | d$pvalues > 1))
            fail("pvalues must be above 0 and at most 1")
    }
    for (v in intersect(c("beta", "varbeta"), read))
        if (!is.numeric(d[[v]]))
            fail(v, " must be numeric")
    if (reads("varbeta") && any(d$varbeta < 0))
        fail("varbeta must be > 0")
    if (reads("sdY") &&
        (!is.numeric(d$sdY) || length(d$sdY) != 1L || !is.finite(d$sdY) || d$sdY <= 0))
        fail("sdY must be a single finite number above 0")
    if (reads("z") &&
        !(is.numeric(d$z) && length(d$z) == length(d$snp) && all(is.finite(d$z))))
        fail("z must be a finite number for each snp")

    ## a warning when no p-value is below warn.minp; a coloc_dataset had it
    ## when it was made. The p-values come from beta and varbeta, or else
    ## pvalues, or else z, which is read here only when it holds numbers
    if (!inherits(d, "coloc_dataset")) {
        p <- if (reads("beta") && reads("varbeta")) pnorm(-abs(d$beta / sqrt(d$varbeta))) * 2
             else if (reads("pvalues")) d$pvalues
             else if ("z" %in% nd && is.numeric(d$z) && !anyNA(d$z)) pnorm(-abs(d$z)) * 2
        if (length(p) && min(p) > warn.minp)
            warning(warningCondition(
                paste0(prefix, "the smallest p-value is ", format.pval(min(p)),
                       ", so this dataset may have no association in the region.",
                       "\nIf this is what you expected, this is not a problem.\nIf this is not as small as you expected, please check you supplied var(beta) and not sd(beta) for the varbeta argument. If that's not the explanation, please check the 02_data vignette."),
                class = "coloc_min_p", call = NULL))
    }
    NULL
}

## for messages: "pvalues is missing", "pvalues and MAF are missing"
.is_missing <- function(v)
    paste(if (length(v) < 2L) v else paste(paste(v[-length(v)], collapse = ", "), "and", v[length(v)]),
          ngettext(length(v), "is missing", "are missing"))

## What each engine reads from a dataset d with the elements nd. An entry
## stops, through fail(), when an element its engine needs is absent, and
## returns the elements whose values check_dataset() then checks. A new engine
## adds an entry here, and the rules for each element's values stay as they
## are. The first entry is check_dataset()'s default.
.engine_reads <- list(
    ## no engine: every element given, apart from z
    none = function(d, nd, fail) setdiff(nd, "z"),
    ## type, then beta and varbeta, or else pvalues, MAF and N
    wakefield = function(d, nd, fail) {
        if (!("type" %in% nd))
            fail('type is missing; set type to "quant" for a quantitative trait or "cc" for a case-control trait')
        ## the type decides what else is read. check_dataset() tests it
        ## afterwards, so a type that is not "cc" is read here as "quant"
        cc <- isTRUE(d$type == "cc")
        if (all(c("beta", "varbeta") %in% nd)) {
            ## a quantitative trait's prior is in units of sdY, which is
            ## estimated from varbeta, MAF and N when not given
            if (cc)
                return(c("type", "beta", "varbeta"))
            if ("sdY" %in% nd)
                return(c("type", "beta", "varbeta", "sdY"))
            if (!all(c("MAF", "N") %in% nd))
                fail("must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
            return(c("type", "beta", "varbeta", "MAF", "N"))
        }
        if (!all(c("pvalues", "MAF") %in% nd))
            fail(.is_missing(setdiff(c("pvalues", "MAF"), nd)),
                 "; without beta and varbeta, coloc needs pvalues, MAF and N")
        if (cc && !("s" %in% nd))
            fail("s is missing; a case-control dataset without beta and varbeta needs s, ",
                 "the proportion of samples who are cases")
        if (!("N" %in% nd))
            fail("N is missing; without beta and varbeta, coloc needs N, the sample size")
        c("type", "pvalues", "MAF", "N", if (cc) "s")
    },
    ## LD and N, then z, or else beta and varbeta
    susie = function(d, nd, fail) {
        if (!("LD" %in% nd))
            fail("LD is missing; SuSiE needs LD, the correlation matrix of the snps")
        if (!("N" %in% nd))
            fail("N is missing; SuSiE needs N, the sample size")
        if ("z" %in% nd)
            return(c("LD", "N", "z"))
        if (!all(c("beta", "varbeta") %in% nd))
            fail(.is_missing(setdiff(c("beta", "varbeta"), nd)),
                 "; without z, SuSiE needs beta and varbeta")
        c("LD", "N", "beta", "varbeta")
    })


## The refusal of anything but a list of dataset elements, for check_dataset()
## and new_coloc_dataset(). With a label, such as 1, the message begins
## "dataset 1: ".
.stop_not_a_dataset <- function(d, dataset = NULL)
    stop(if (is.null(dataset)) "the dataset" else paste0("dataset ", dataset, ": this"),
         " is ", .what_is(d), ", not a list of coloc_dataset elements (snp, beta, varbeta, ...)",
         call. = FALSE)

## TRUE when LD is a correlation matrix, to within tol: no missing values,
## values in [-1, 1], 1 on the diagonal, and symmetric. The default passes an
## LD worked out in 64-bit arithmetic or stored as 32-bit floats, but one
## worked out in 32-bit arithmetic can be out by about 2e-5, and is refused.
## The transpose is compared in blocks of 500 columns, so a large LD isn't
## copied whole.
.is_correlation <- function(LD, tol = 1e-6) {
    if (anyNA(LD))
        return(FALSE)
    r <- range(LD)
    if (r[1] < -1 - tol || r[2] > 1 + tol)
        return(FALSE)
    if (any(abs(diag(LD) - 1) > tol))
        return(FALSE)
    for (j in split(seq_len(ncol(LD)), ceiling(seq_len(ncol(LD)) / 500)))
        if (max(abs(LD[, j, drop = FALSE] - t(LD[j, , drop = FALSE]))) > tol)
            return(FALSE)
    TRUE
}

## What x is, for a message: "a numeric vector", "a data.frame", "a SuSiE fit"
.what_is <- function(x) {
    if (is.null(x))
        return("NULL")
    what <- if (inherits(x, "susie")) "SuSiE fit"
            else if (is.object(x) || !is.atomic(x)) class(x)[1]
            else paste(mode(x), if (is.matrix(x)) "matrix" else "vector")
    paste(if (grepl("^[aeiou]", what)) "an" else "a", what)
}
