## log Bayes factors for two signals over five snps
lbf <- matrix(c(0.1, 3.2, 0.4, -0.5, 1.0,
                2.5, 0.3, 0.0, 6.1, -1.2), nrow = 2, byrow = TRUE,
              dimnames = list(NULL, paste0("s", 1:5)))

## a coloc_finemap as the validator sees it, without the constructor's help
raw_finemap <- function(...) {
    x <- list(bf = lbf, idx = 1:2, snps = data.frame(snp = colnames(lbf)),
              engine = "external", args = list())
    x[names(list(...))] <- list(...)
    structure(x, class = "coloc_finemap")
}
## the message validate_coloc_finemap() stops with
refusal <- function(x) conditionMessage(expect_error(validate_coloc_finemap(x)))

test_that("as_coloc_finemap() makes a coloc_finemap from a matrix of log Bayes factors", {
    f <- as_coloc_finemap(lbf)
    expect_s3_class(f, "coloc_finemap")
    expect_identical(f$bf, lbf)
    expect_identical(f$idx, 1:2)
    expect_identical(f$snps, data.frame(snp = paste0("s", 1:5)))
    expect_identical(f$engine, "external")
    expect_identical(f$args, list())
})

test_that("an integer matrix is stored as double", {
    m <- matrix(1:10, nrow = 2, dimnames = list(NULL, paste0("s", 1:5)))
    f <- as_coloc_finemap(m)
    expect_identical(storage.mode(f$bf), "double")
    expect_equal(f$bf, m)
})

test_that("bf has no row names, as idx names its rows", {
    m <- lbf
    rownames(m) <- c("L1", "L2")
    expect_null(rownames(as_coloc_finemap(m)$bf))
    f <- new_coloc_finemap(m, c(3L, 7L), data.frame(snp = colnames(lbf)), "susie")
    expect_null(rownames(f$bf))
    expect_identical(unname(f$bf), unname(lbf))
})

test_that("as_coloc_finemap() makes a one-row coloc_finemap from a named vector", {
    v <- c(s1 = 0.1, s2 = 3.2, s3 = 0.4)
    f <- as_coloc_finemap(v)
    expect_identical(f$bf, matrix(v, nrow = 1, dimnames = list(NULL, names(v))))
    expect_identical(f$idx, 1L)
    expect_identical(f$snps, data.frame(snp = names(v)))
    expect_identical(as_coloc_finemap(c(s1 = 1L, s2 = 2L))$bf,
                     matrix(c(1, 2), nrow = 1, dimnames = list(NULL, c("s1", "s2"))))
})

test_that("a coloc_finemap is returned as it is, even one changed since it was made", {
    f <- as_coloc_finemap(lbf)
    expect_identical(as_coloc_finemap(f), f)
    f$idx <- c(4, 9)                    # still a coloc_finemap: $<- keeps the class
    expect_identical(as_coloc_finemap(f), f)
    broken <- f
    broken$idx <- 1                     # one idx for two rows
    expect_identical(as_coloc_finemap(broken), broken)
})

test_that("as_coloc_finemap() refuses a matrix or vector without snp ids", {
    e <- expect_error(as_coloc_finemap(unname(lbf)))
    expect_identical(conditionMessage(e), "the matrix has no column names; name its columns by snp id, or give snps")
    expect_null(conditionCall(e))
    e <- expect_error(as_coloc_finemap(c(0.1, 3.2)))
    expect_identical(conditionMessage(e), "the vector has no names; name it by snp id, or give snps")
})

test_that("as_coloc_finemap() says what it takes when given anything else", {
    takes <- "; as_coloc_finemap() takes a matrix of log Bayes factors, a numeric vector or a coloc_finemap"
    e <- expect_error(as_coloc_finemap(list(bf = lbf)))
    expect_identical(conditionMessage(e), paste0("x is a list", takes))
    expect_null(conditionCall(e))
    e <- expect_error(as_coloc_finemap(as.data.frame(lbf)))
    expect_identical(conditionMessage(e), paste0("x is a data.frame", takes))
    e <- expect_error(as_coloc_finemap(c(s1 = TRUE)))
    expect_identical(conditionMessage(e), paste0("x is a logical vector", takes))
})

test_that("a coloc_finemap with no signal passes", {
    none <- lbf[0, , drop = FALSE]
    f <- as_coloc_finemap(none)
    expect_identical(dim(f$bf), c(0L, 5L))
    expect_identical(f$idx, integer(0))
    expect_identical(f$snps$snp, colnames(lbf))
    f <- new_coloc_finemap(none, integer(0), data.frame(snp = colnames(lbf)), "susie")
    expect_identical(f$idx, integer(0))
})

test_that("a null column stays in bf, last, and is left out of snps", {
    m <- cbind(null = c(0, 0), lbf)
    f <- as_coloc_finemap(m)
    expect_identical(colnames(f$bf), c(paste0("s", 1:5), "null"))
    expect_identical(f$bf[, "null"], c(0, 0))
    expect_identical(f$bf[, 1:5], lbf)
    expect_identical(f$snps$snp, paste0("s", 1:5))
    ## in the middle, with the snps kept in their order
    m <- cbind(lbf[, 1:2], null = c(1, 2), lbf[, 3:5])
    f <- as_coloc_finemap(m)
    expect_identical(colnames(f$bf), c(paste0("s", 1:5), "null"))
    expect_identical(f$bf[, "null"], c(1, 2))
})

## lbf's snp ids, and lbf without them
ids <- paste0("s", 1:5)
bare <- unname(lbf)

test_that("as_coloc_finemap() takes the snp ids of a matrix without column names from a snps data.frame", {
    f <- as_coloc_finemap(bare, snps = data.frame(snp = ids))
    expect_identical(f$bf, lbf)
    expect_identical(f$snps, data.frame(snp = ids))
})

test_that("as_coloc_finemap() takes snps as a character vector of snp ids", {
    f <- as_coloc_finemap(bare, snps = ids)
    expect_identical(f$bf, lbf)
    expect_identical(f$snps, data.frame(snp = ids))
})

test_that("as_coloc_finemap() keeps the other columns of a snps data.frame", {
    snps <- data.frame(snp = ids, position = c(100L, 250L, 300L, 410L, 555L))
    expect_identical(as_coloc_finemap(bare, snps = snps)$snps, snps)
    expect_identical(as_coloc_finemap(lbf, snps = snps)$snps, snps)
})

test_that("as_coloc_finemap() takes snps for a vector, named or not", {
    v <- c(s1 = 0.1, s2 = 3.2, s3 = 0.4)
    f <- as_coloc_finemap(v, snps = names(v))
    expect_identical(f$bf, matrix(v, nrow = 1, dimnames = list(NULL, names(v))))
    expect_identical(as_coloc_finemap(unname(v), snps = names(v)), f)
    expect_identical(as_coloc_finemap(unname(v), snps = data.frame(snp = names(v))), f)
})

test_that("as_coloc_finemap() refuses snps with the wrong number of snps", {
    e <- expect_error(as_coloc_finemap(bare, snps = ids[1:4]))
    expect_identical(conditionMessage(e), "snps has 4 snps but x has 5 log Bayes factors per signal")
    expect_null(conditionCall(e))
    e <- expect_error(as_coloc_finemap(c(0.1, 3.2), snps = data.frame(snp = ids[1:3])))
    expect_identical(conditionMessage(e), "snps has 3 snps but x has 2 log Bayes factors per signal")
})

test_that("as_coloc_finemap() refuses snps that differ from the snp ids of x", {
    differ <- "snps$snp must be the column names of bf, in the same order, without null"
    e <- expect_error(as_coloc_finemap(lbf, snps = c("s1", "s2", "s3", "s4", "s9")))
    expect_identical(conditionMessage(e), differ)
    e <- expect_error(as_coloc_finemap(lbf, snps = rev(ids)))
    expect_identical(conditionMessage(e), differ)
    e <- expect_error(as_coloc_finemap(lbf, snps = ids[1:4]))
    expect_identical(conditionMessage(e), differ)
    ## a null column stays out of snps
    f <- as_coloc_finemap(cbind(lbf, null = c(0, 0)), snps = ids)
    expect_identical(f$snps, data.frame(snp = ids))
})

test_that("as_coloc_finemap() refuses snps that are neither snp ids nor a data.frame with a snp column", {
    takes <- "snps must be a character vector of snp ids or a data.frame with a snp column; it is "
    e <- expect_error(as_coloc_finemap(bare, snps = factor(ids)))
    expect_identical(conditionMessage(e), paste0(takes, "a factor"))
    expect_null(conditionCall(e))
    e <- expect_error(as_coloc_finemap(bare, snps = 1:5))
    expect_identical(conditionMessage(e), paste0(takes, "a numeric vector"))
    e <- expect_error(as_coloc_finemap(bare, snps = as.list(ids)))
    expect_identical(conditionMessage(e), paste0(takes, "a list"))
    e <- expect_error(as_coloc_finemap(lbf, snps = factor(ids)))
    expect_identical(conditionMessage(e), paste0(takes, "a factor"))
    e <- expect_error(as_coloc_finemap(bare, snps = data.frame(id = ids)))
    expect_identical(conditionMessage(e), "snps must have a snp column")
})

test_that("as_coloc_finemap() takes idx, one per signal", {
    expect_identical(as_coloc_finemap(lbf, idx = c(3, 7))$idx, c(3L, 7L))
    expect_identical(as_coloc_finemap(c(s1 = 0.1, s2 = 3.2), idx = 4)$idx, 4L)
    expect_identical(as_coloc_finemap(bare, snps = ids, idx = 2:1)$idx, 2:1)
})

test_that("as_coloc_finemap() refuses an idx that is not one whole number above 0 per signal, with no repeat", {
    e <- expect_error(as_coloc_finemap(lbf, idx = 1:3))
    expect_identical(conditionMessage(e), "idx has 3 values but bf has 2 rows")
    whole <- "idx must be whole numbers above 0, with no missing value and no repeat"
    e <- expect_error(as_coloc_finemap(lbf, idx = c(1, 1)))
    expect_identical(conditionMessage(e), whole)
    e <- expect_error(as_coloc_finemap(lbf, idx = c(0, 1)))
    expect_identical(conditionMessage(e), whole)
})

test_that("as_coloc_finemap() refuses snps or idx for a coloc_finemap", {
    f <- as_coloc_finemap(lbf)
    already <- "x is already a coloc_finemap; snps and idx apply only to a matrix or a vector"
    e <- expect_error(as_coloc_finemap(f, snps = ids))
    expect_identical(conditionMessage(e), already)
    expect_null(conditionCall(e))
    e <- expect_error(as_coloc_finemap(f, idx = 1:2))
    expect_identical(conditionMessage(e), already)
})

test_that("as_coloc_finemap() turns a factor snp column into character, as R before 4.0 makes one by default", {
    snps <- data.frame(snp = factor(ids), position = 1:5)
    f <- as_coloc_finemap(bare, snps = snps)
    expect_identical(f$snps, data.frame(snp = ids, position = 1:5))
    expect_identical(colnames(f$bf), ids)
    expect_identical(as_coloc_finemap(lbf, snps = snps)$snps, f$snps)
})

test_that("as_coloc_finemap() refuses an argument it does not take", {
    v <- c(s1 = 0.1, s2 = 3.2)
    e <- expect_error(as_coloc_finemap(lbf, idxx = c(5, 9)))
    expect_identical(conditionMessage(e), "unused argument: idxx")
    expect_null(conditionCall(e))
    e <- expect_error(as_coloc_finemap(v, idxx = 4, foo = 1))
    expect_identical(conditionMessage(e), "unused arguments: idxx, foo")
    e <- expect_error(as_coloc_finemap(as_coloc_finemap(lbf), snp = ids))
    expect_identical(conditionMessage(e), "unused argument: snp")
})

test_that("the validator refuses a coloc_finemap that is not a list", {
    expect_identical(refusal(structure(1:3, class = "coloc_finemap")),
                     "coloc_finemap must be a list; it is a numeric vector")
})

test_that("the validator reads each element by its exact name", {
    ## $ would take bfx for bf, idx2 for idx, and so on
    msgs <- c(bf = "bf must be a numeric matrix of log Bayes factors, signals x snps; it is NULL",
              idx = "idx must be whole numbers above 0, with no missing value and no repeat",
              snps = "snps must be a data.frame; it is NULL",
              engine = "engine must be one of: wakefield, susie, external",
              args = "args must be a list; it is NULL")
    for (nm in names(msgs)) {
        x <- raw_finemap()
        names(x)[names(x) == nm] <- paste0(nm, "x")
        expect_identical(refusal(x), msgs[[nm]])
    }
})

test_that("the validator refuses a bf that is not a numeric matrix with snp ids as column names", {
    expect_identical(refusal(raw_finemap(bf = as.data.frame(lbf))),
                     "bf must be a numeric matrix of log Bayes factors, signals x snps; it is a data.frame")
    expect_identical(conditionMessage(expect_error(as_coloc_finemap(matrix("a", 1, 1, dimnames = list(NULL, "s1"))))),
                     "bf must be a numeric matrix of log Bayes factors, signals x snps; it is a character matrix")
    expect_identical(refusal(raw_finemap(bf = unname(lbf))), "bf must have snp ids as column names")
    m <- lbf
    colnames(m)[3] <- ""
    expect_identical(conditionMessage(expect_error(as_coloc_finemap(m))),
                     "bf has an empty column name; name every column by snp id")
    colnames(m)[3] <- NA
    expect_identical(conditionMessage(expect_error(as_coloc_finemap(m))),
                     "bf has an empty column name; name every column by snp id")
    expect_identical(conditionMessage(expect_error(as_coloc_finemap(c(s1 = 1, 2)))),
                     "bf has an empty column name; name every column by snp id")
})

test_that("the validator refuses a bf without snps", {
    e <- expect_error(as_coloc_finemap(matrix(numeric(0), 1, 0)))
    expect_identical(conditionMessage(e), "bf has no snps")
    e <- expect_error(as_coloc_finemap(c(null = 0)))
    expect_identical(conditionMessage(e), "bf has no snps")
})

test_that("the validator refuses a missing or infinite log Bayes factor", {
    for (v in c(NA, NaN, Inf, -Inf)) {
        m <- lbf
        m[2, 4] <- v
        e <- expect_error(as_coloc_finemap(m))
        expect_identical(conditionMessage(e), "bf contains missing or infinite values")
    }
})

test_that("the validator refuses a snp id twice, which also refuses two null columns", {
    m <- lbf
    colnames(m)[5] <- "s2"
    e <- expect_error(as_coloc_finemap(m))
    expect_identical(conditionMessage(e), "snp id s2 is in bf twice")
    e <- expect_error(as_coloc_finemap(cbind(lbf, null = 0, null = 0)))
    expect_identical(conditionMessage(e), "snp id null is in bf twice")
})

test_that("the validator refuses a null column that is not last", {
    x <- raw_finemap(bf = cbind(null = 0, lbf))
    expect_identical(refusal(x), "the null column of bf must be the last")
})

test_that("the validator refuses an idx that is not one whole number above 0 per signal, with no repeat", {
    msg <- "idx must be whole numbers above 0, with no missing value and no repeat"
    for (idx in list(c(1, 1), c(0, 1), c(-1, 2), c(1.5, 2), c(1, NA), c(1, Inf), c(1, 2^31),
                     c("1", "2")))
        expect_identical(refusal(raw_finemap(idx = idx)), msg)
    expect_identical(refusal(raw_finemap(idx = 1L)), "idx has 1 value but bf has 2 rows")
    ## whole numbers stored as double pass, and are stored as integer
    f <- new_coloc_finemap(lbf, c(3, 7), data.frame(snp = colnames(lbf)), "susie")
    expect_identical(f$idx, c(3L, 7L))
})

test_that("the validator refuses snps that are not a data.frame of the snp ids in bf's order", {
    expect_identical(refusal(raw_finemap(snps = list(snp = colnames(lbf)))),
                     "snps must be a data.frame; it is a list")
    expect_identical(refusal(raw_finemap(snps = data.frame(id = colnames(lbf)))),
                     "snps must have a snp column")
    msg <- "snps$snp must be the column names of bf, in the same order, without null"
    expect_identical(refusal(raw_finemap(snps = data.frame(snp = rev(colnames(lbf))))), msg)
    expect_identical(refusal(raw_finemap(snps = data.frame(snp = colnames(lbf)[1:4]))), msg)
    expect_identical(refusal(raw_finemap(snps = data.frame(snp = factor(colnames(lbf))))),
                     "snps$snp must be a character vector of snp ids; it is a factor")
    expect_identical(refusal(raw_finemap(bf = cbind(lbf, null = 0),
                                         snps = data.frame(snp = c(colnames(lbf), "null")))), msg)
    ## more columns are fine
    s <- data.frame(snp = colnames(lbf), position = 1:5 * 100L, z = rnorm(5))
    expect_identical(new_coloc_finemap(lbf, 1:2, s, "wakefield")$snps, s)
})

test_that("the validator refuses an unknown engine and args that are not a list", {
    expect_identical(refusal(raw_finemap(engine = "abf")),
                     "engine must be one of: wakefield, susie, external")
    expect_identical(refusal(raw_finemap(engine = c("susie", "wakefield"))),
                     "engine must be one of: wakefield, susie, external")
    expect_identical(refusal(raw_finemap(args = 1)), "args must be a list; it is a numeric vector")
    expect_null(conditionCall(expect_error(validate_coloc_finemap(raw_finemap(args = 1)))))
})

test_that("print() shows the engine, the counts and each signal's top snp", {
    expect_identical(capture.output(print(as_coloc_finemap(lbf))),
                     c("coloc finemap: external, 2 signals, 5 snps",
                       "  signal 1 (idx 1): top snp s2, log BF 3.2",
                       "  signal 2 (idx 2): top snp s4, log BF 6.1"))
    f <- new_coloc_finemap(lbf[2, , drop = FALSE] * 1.23456, 4L, data.frame(snp = colnames(lbf)), "susie")
    expect_identical(capture.output(print(f)),
                     c("coloc finemap: susie, 1 signal, 5 snps",
                       "  signal 1 (idx 4): top snp s4, log BF 7.53"))
    expect_identical(capture.output(print(as_coloc_finemap(lbf[0, , drop = FALSE]))),
                     c("coloc finemap: external, 0 signals, 5 snps",
                       "  no signals"))
    expect_identical(capture.output(print(as_coloc_finemap(c(s1 = 2)))),
                     c("coloc finemap: external, 1 signal, 1 snp",
                       "  signal 1 (idx 1): top snp s1, log BF 2"))
    capture.output(expect_invisible(print(as_coloc_finemap(lbf))))
})

test_that("print() leaves the null column out of the top snp and the count of snps", {
    f <- as_coloc_finemap(cbind(lbf, null = c(9, 9)))
    expect_identical(capture.output(print(f)),
                     c("coloc finemap: external, 2 signals, 5 snps",
                       "  signal 1 (idx 1): top snp s2, log BF 3.2",
                       "  signal 2 (idx 2): top snp s4, log BF 6.1"))
})
