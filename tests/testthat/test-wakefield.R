D1 <- coloc_test_data$D1

## Wakefield's log ABF from z, the variance of the effect estimate V and the
## prior standard deviation of the effect
labf <- function(z, V, sd_prior) {
    r <- sd_prior^2 / (sd_prior^2 + V)
    0.5 * (log(1 - r) + r * z^2)
}

test_that("with beta, varbeta and sdY, the log ABF follows Wakefield's formula", {
    got <- wakefield_bf(as_coloc_dataset(D1))
    z <- D1$beta / sqrt(D1$varbeta)
    expect_equal(got$z, z)
    expect_equal(got$V, D1$varbeta)
    expect_equal(got$r, (0.15 * D1$sdY)^2 / ((0.15 * D1$sdY)^2 + D1$varbeta))
    expect_equal(got$lABF, labf(z, D1$varbeta, 0.15 * D1$sdY))
})

test_that("a case-control trait with beta and varbeta has its prior on the log odds ratio scale", {
    got <- wakefield_bf(as_coloc_dataset(as_cc(D1)))
    expect_equal(got$lABF, labf(D1$beta / sqrt(D1$varbeta), D1$varbeta, 0.2))
})

test_that("without sdY, wakefield_bf() estimates it from varbeta, MAF and N, and warns", {
    d <- D1
    d$sdY <- NULL
    d <- as_coloc_dataset(d)
    expect_warning(got <- wakefield_bf(d),
                   "estimating sdY from maf and varbeta, please directly supply sdY if known",
                   fixed = TRUE)
    ## var(beta-hat) = var(Y) / (2 N maf (1 - maf)), fitted through the origin
    x <- 1 / D1$varbeta
    y <- 2 * D1$N * D1$MAF * (1 - D1$MAF)
    sdY <- sqrt(sum(x * y) / sum(x^2))
    expect_equal(got$lABF, labf(D1$beta / sqrt(D1$varbeta), D1$varbeta, 0.15 * sdY))
})

test_that("from p-values, z and V are approximated from pvalues, MAF and N", {
    d <- to_pvalues(D1)
    got <- wakefield_bf(as_coloc_dataset(d))
    z <- qnorm(d$pvalues / 2, lower.tail = FALSE)
    V <- 1 / (2 * d$N * d$MAF * (1 - d$MAF))
    expect_equal(got$z, z)
    expect_equal(got$V, V)
    expect_equal(got$lABF, labf(z, V, 0.15))
})

test_that("from p-values, a case-control trait's V also takes s", {
    d <- to_pvalues(as_cc(D1))
    got <- wakefield_bf(as_coloc_dataset(d))
    V <- 1 / (2 * d$N * d$MAF * (1 - d$MAF) * d$s * (1 - d$s))
    expect_equal(got$V, V)
    expect_equal(got$lABF, labf(qnorm(d$pvalues / 2, lower.tail = FALSE), V, 0.2))
})

test_that("wakefield_bf() returns a data.frame with one row per snp", {
    got <- wakefield_bf(as_coloc_dataset(D1))
    expect_identical(class(got), "data.frame")
    expect_named(got, c("snp", "position", "z", "V", "r", "lABF"))
    expect_identical(got$snp, D1$snp)
    expect_identical(got$position, D1$position)
    ## the row names are row numbers, whatever the dataset's vectors are named
    named <- D1
    names(named$beta) <- D1$snp
    names(named$position) <- D1$snp
    expect_identical(rownames(wakefield_bf(as_coloc_dataset(named))),
                     as.character(seq_along(D1$snp)))
})

test_that("position is NA for a dataset without one", {
    for (d in list(no_position(D1), no_position(to_pvalues(D1)))) {
        got <- wakefield_bf(as_coloc_dataset(d))
        expect_identical(got$position, rep(NA_integer_, length(D1$snp)))
    }
    ## no snps and no position: no rows, where a single NA position would fail
    none <- no_position(D1)
    none$LD <- NULL
    for (v in c("beta", "varbeta", "MAF", "snp"))
        none[[v]] <- none[[v]][0]
    none <- suppressWarnings(as_coloc_dataset(none))
    expect_identical(nrow(wakefield_bf(none)), 0L)
})

test_that("effect_priors sets the prior for the dataset's type of trait, on both routes", {
    wide <- c(quant = 0.3, cc = 0.4)
    for (d in list(D1, as_cc(D1), to_pvalues(D1), to_pvalues(as_cc(D1)))) {
        d <- as_coloc_dataset(d)
        default <- wakefield_bf(d)
        sdY <- if (is.null(d$sdY)) 1 else d$sdY
        expect_equal(wakefield_bf(d, effect_priors = wide)$lABF,
                     labf(default$z, default$V, wide[[d$type]] * sdY))
        ## the prior for the other type of trait changes nothing
        other <- if (d$type == "quant") c(quant = 0.15, cc = 0.4) else c(quant = 0.3, cc = 0.2)
        expect_identical(wakefield_bf(d, effect_priors = other), default)
    }
})

test_that("wakefield_bf() refuses effect_priors without a positive number for each type", {
    D <- as_coloc_dataset(D1)
    for (bad in list(c(quant = 0.15), c(0.15, 0.2), c(quant = 0, cc = 0.2),
                     c(quant = "0.15", cc = "0.2"))) {
        e <- expect_error(wakefield_bf(D, effect_priors = bad),
                          "effect_priors should be positive numbers named quant and cc", fixed = TRUE)
        expect_null(conditionCall(e))
    }
})

test_that("wakefield_bf() checks a coloc_dataset again, as it may have changed since it was made", {
    changed <- as_coloc_dataset(D1)
    changed$beta[1] <- NA
    e <- expect_error(wakefield_bf(changed))
    expect_identical(conditionMessage(e), "beta contains missing values")
    expect_null(conditionCall(e))
})

test_that("wakefield_bf() refuses a coloc_dataset without what Wakefield needs", {
    ## made with only what SuSiE needs
    z_only <- as_coloc_dataset(list(snp = D1$snp, z = D1$beta / sqrt(D1$varbeta),
                                    LD = D1$LD, N = D1$N))
    e <- expect_error(wakefield_bf(z_only))
    expect_identical(conditionMessage(e),
                     'type is missing; set type to "quant" for a quantitative trait or "cc" for a case-control trait')
    no_n <- to_pvalues(D1)
    no_n$N <- NULL
    e <- expect_error(wakefield_bf(as_coloc_dataset(no_n)))
    expect_identical(conditionMessage(e),
                     "N is missing; without beta and varbeta, coloc needs pvalues, MAF and N, the sample size")
})

test_that("a weak dataset gets the smallest p-value warning when it is made, and not from wakefield_bf()", {
    weak <- D1
    weak$beta <- weak$beta / 100
    expect_warning(weak <- as_coloc_dataset(weak), class = "coloc_min_p")
    expect_no_warning(wakefield_bf(weak))
})

test_that("Var.data() and Var.data.cc() give the variance of the effect estimate", {
    expect_equal(Var.data(0.2, 1000), 1 / (2 * 1000 * 0.2 * 0.8))
    expect_equal(Var.data.cc(0.2, 1000, 0.3), 1 / (2 * 1000 * 0.2 * 0.8 * 0.3 * 0.7))
})

test_that("sdY.est() recovers a trait's standard deviation, and warns that it is an estimate", {
    set.seed(42)
    X <- sapply(runif(100, 0.05, 0.5), function(maf) rbinom(400, 2, maf))
    Y <- rnorm(nrow(X), mean = X[, 1], sd = 4) + rnorm(nrow(X), sd = 2)
    vbeta <- apply(X, 2, function(x) vcov(lm(Y ~ x))[2, 2])
    expect_warning(est <- sdY.est(vbeta = vbeta, maf = colMeans(X) / 2, n = nrow(X)),
                   "estimating sdY from maf and varbeta", fixed = TRUE)
    expect_lt(abs(est - sd(Y)), 0.1)
})

test_that("sdY.est() stops when its estimate is negative", {
    expect_error(suppressWarnings(sdY.est(vbeta = c(-1, -2), maf = c(0.2, 0.3), n = 100)),
                 "estimated sdY is negative", fixed = TRUE)
})
