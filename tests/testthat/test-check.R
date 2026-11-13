D1 <- coloc_test_data$D1

B1=D1
C1=D1
C1[c("type","s")]=list(type="cc",s=0.5)
B1[c("type","s")]=list(type="cc",s=2)

test_that("missing or bad required elements throws error", {
  expect_null(check_dataset(D1, req = names(D1)))
  e <- expect_error(check_dataset(D1, req = "test"))
  expect_identical(conditionMessage(e), "missing required element(s) test")
  ## type="cc" and s
  expect_null(check_dataset(C1))
  e <- expect_error(check_dataset(B1))
  expect_identical(conditionMessage(e), "s must be between 0 and 1")
  e <- expect_error(check_dataset(list(), ""))
  expect_identical(conditionMessage(e), "missing required element(s) snp")
  e <- expect_error(check_dataset(list(beta=1,p=2,type="blah"), ""))
  expect_identical(conditionMessage(e), "missing required element(s) snp")
})

test_that("issue 79", {
  d2=list(snp=letters[1:5],
          position=1:5,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          type="quant",
          MAF=rep(0.5,5))
  e <- expect_error(check_dataset(d2, 2, engine = "wakefield"))
  expect_identical(conditionMessage(e),
                   "dataset 2: must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
})

test_that("issue 160", {
  d2=list(snp=letters[1:5],
          position=1:5,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          N=rep(150000,5),
          type="quant",
          MAF=rep(0.5,5))
  e <- expect_error(check_dataset(d2, 2))
  expect_identical(conditionMessage(e),
                   "dataset 2: N must be a single number; coloc does not take per-snp sample sizes")
})

test_that("Infinite values in beta/varbeta triggers specific error", {
  D1_beta_inf <- D1
  D1_beta_inf$beta[[1]] <- Inf
  e <- expect_error(check_dataset(D1_beta_inf))
  expect_identical(conditionMessage(e), "beta has infinite values")

  D1_varbeta_inf <- D1
  D1_varbeta_inf$varbeta[[1]] <- Inf
  e <- expect_error(check_dataset(D1_varbeta_inf))
  expect_identical(conditionMessage(e), "varbeta has infinite values")

  D1_varbeta_inf <- D1
  D1_varbeta_inf$varbeta[[1]] <- 0
  e <- expect_error(check_dataset(D1_varbeta_inf))
  expect_identical(conditionMessage(e), "varbeta must be > 0")
})

## ---- check_dataset(), one test per check -----------------------------------
## Most of check_dataset()'s checks have a case here: a dataset that fails it,
## which check_dataset() refuses. Other cases
## are valid datasets, datasets that fail two checks, and the arguments engine,
## req and warn.minp. Each case runs with the label 1, as coloc() checks a dataset, and
## with none, as a dataset is made. A refused case holds the message
## check_dataset() stops with, as coloc() shows it with the label 1, so the
## test shows which check refused the dataset. A case that gives a warning holds
## its first line, which says what is wrong.

## the datasets the cases are made from: D1 as given (beta and varbeta), D1 as
## p-values, MAF and N, and D1 as z-scores, with only what SuSiE reads
pv1 <- to_pvalues(D1)
z1 <- list(snp = D1$snp, z = D1$beta / sqrt(D1$varbeta), LD = D1$LD, N = D1$N)
wakefield <- list(engine = "wakefield")
susie <- list(engine = "susie")
## each case: the dataset, any arguments besides the dataset and its suffix,
## the message check_dataset() refuses it with, and the first line of the
## warning it gives
cases <- list(
    ## valid datasets
    "a valid dataset passes" =
        list(d = D1),
    "a dataset given as p-values passes" =
        list(d = pv1),
    "a case-control trait with beta and varbeta needs neither s nor sdY" =
        list(d = local({
                 d <- D1
                 d$type <- "cc"
                 d$sdY <- NULL
                 d
             })),
    "a dataset without LD passes" =
        list(d = local({
                 d <- D1
                 d$LD <- NULL
                 d
             })),
    "an element check_dataset() doesn't know may hold missing values" =
        list(d = local({
                 d <- D1
                 d$extra <- NA
                 d
             })),
    ## a dataset for each check
    "the dataset is a list" =
        list(d = 1:3,
             error = "dataset 1: this is a numeric vector, not a list of coloc_dataset elements (snp, beta, varbeta, ...)"),
    "the required elements are there" =
        list(d = local({
                 d <- D1
                 d$snp <- NULL
                 d
             }),
             error = "dataset 1: missing required element(s) snp"),
    "snp is the only element every dataset needs" =
        list(d = list(snp = D1$snp)),
    "with no engine, a dataset needs no type" =
        list(d = local({
                 d <- D1
                 d$type <- NULL
                 d
             })),
    "no missing values in beta" =
        list(d = local({
                 d <- D1
                 d$beta[1] <- NA
                 d
             }),
             error = "dataset 1: beta contains missing values"),
    "no missing values in MAF" =
        list(d = local({
                 d <- D1
                 d$MAF[1] <- NA
                 d
             }),
             error = "dataset 1: MAF contains missing values"),
    "type is quant or cc" =
        list(d = local({
                 d <- D1
                 d$type <- "binary"
                 d
             }),
             error = "dataset 1: type must be quant or cc"),
    "the snp ids are unique" =
        list(d = local({
                 d <- D1
                 d$snp[2] <- d$snp[1]
                 d
             }),
             error = "dataset 1: duplicated snps found"),
    "snp is not a factor" =
        list(d = local({
                 d <- D1
                 d$snp <- factor(d$snp)
                 d
             }),
             error = "dataset 1: snp should be a character vector but is a factor"),
    "MAF is above 0" =
        list(d = local({
                 d <- D1
                 d$MAF[1] <- 0
                 d
             }),
             error = "dataset 1: MAF should be a numeric, strictly >0 & <1"),
    "MAF is below 1" =
        list(d = local({
                 d <- D1
                 d$MAF[1] <- 1
                 d
             }),
             error = "dataset 1: MAF should be a numeric, strictly >0 & <1"),
    "MAF is numeric" =
        list(d = local({
                 d <- D1
                 d$MAF <- as.character(d$MAF)
                 d
             }),
             error = "dataset 1: MAF should be a numeric, strictly >0 & <1"),
    "the vectors have the same length" =
        list(d = local({
                 d <- D1
                 d$beta <- d$beta[-1]
                 d
             }),
             error = "dataset 1: snp has 500 values but beta has 499"),
    "N is a single number" =
        list(d = local({
                 d <- D1
                 d$N <- c(1000, 1000)
                 d
             }),
             error = "dataset 1: N must be a single number; coloc does not take per-snp sample sizes"),
    "N is numeric" =
        list(d = local({
                 d <- D1
                 d$N <- as.character(d$N)
                 d
             }),
             error = "dataset 1: N must be a single number; coloc does not take per-snp sample sizes"),
    "s is between 0 and 1" =
        list(d = local({
                 d <- D1
                 d$type <- "cc"
                 d$s <- 1.5
                 d
             }),
             error = "dataset 1: s must be between 0 and 1"),
    "s is above 0" =
        list(d = local({
                 d <- D1
                 d$type <- "cc"
                 d$s <- 0
                 d
             }),
             error = "dataset 1: s must be between 0 and 1"),
    "s is below 1" =
        list(d = local({
                 d <- D1
                 d$type <- "cc"
                 d$s <- 1
                 d
             }),
             error = "dataset 1: s must be between 0 and 1"),
    "s is numeric" =
        list(d = local({
                 d <- D1
                 d$type <- "cc"
                 d$s <- "0.5"
                 d
             }),
             error = "dataset 1: s must be between 0 and 1"),
    "wakefield: beta without varbeta needs p-values and MAF" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$varbeta <- NULL
                 d
             }),
             error = "dataset 1: pvalues is missing; without beta and varbeta, coloc needs pvalues, MAF and N"),
    "wakefield: without beta and varbeta, there are p-values" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$beta <- NULL
                 d$varbeta <- NULL
                 d
             }),
             error = "dataset 1: pvalues is missing; without beta and varbeta, coloc needs pvalues, MAF and N"),
    "wakefield: without beta and varbeta, there is MAF" =
        list(args = wakefield, d = local({
                 d <- pv1
                 d$MAF <- NULL
                 d
             }),
             error = "dataset 1: MAF is missing; without beta and varbeta, coloc needs pvalues, MAF and N"),
    "wakefield: a case-control trait given as p-values has s" =
        list(args = wakefield, d = local({
                 d <- pv1
                 d$type <- "cc"
                 d
             }),
             error = "dataset 1: s is missing; a case-control dataset without beta and varbeta needs s, the proportion of samples who are cases"),
    "wakefield: a trait given as p-values has N" =
        list(args = wakefield, d = local({
                 d <- pv1
                 d$N <- NULL
                 d
             }),
             error = "dataset 1: N is missing; without beta and varbeta, coloc needs N, the sample size"),
    "a trait given as p-values has N above 1" =
        list(d = local({
                 d <- pv1
                 d$N <- 1
                 d
             }),
             error = "dataset 1: N must be finite and > 1"),
    "beta and varbeta are finite" =
        list(d = local({
                 d <- D1
                 d$beta[1] <- Inf
                 d
             }),
             error = "dataset 1: beta has infinite values"),
    "varbeta has no zeros" =
        list(d = local({
                 d <- D1
                 d$varbeta[1] <- 0
                 d
             }),
             error = "dataset 1: varbeta must be > 0"),
    "varbeta has no negative values" =
        list(d = local({
                 d <- D1
                 d$varbeta[1] <- -1
                 d
             }),
             error = "dataset 1: varbeta must be > 0"),
    "a warning when no p-value is small" =
        list(d = local({
                 d <- D1
                 d$beta <- d$beta / 100
                 d
             }),
             warning = "dataset 1: the smallest p-value is 0.93383, so this dataset may have no association in the region."),
    "a warning when no p-value is small, given as p-values" =
        list(d = local({
                 d <- pv1
                 d$pvalues[] <- 0.5
                 d
             }),
             warning = "dataset 1: the smallest p-value is 0.5, so this dataset may have no association in the region."),
    "wakefield: a quantitative trait without sdY has MAF" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$sdY <- NULL
                 d$MAF <- NULL
                 d
             }),
             error = "dataset 1: must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated"),
    "wakefield: a quantitative trait without sdY has N" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$sdY <- NULL
                 d$N <- NULL
                 d
             }),
             error = "dataset 1: must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated"),
    "LD is square" =
        list(d = local({
                 d <- D1
                 d$LD <- d$LD[, -1]
                 d
             }),
             error = "dataset 1: LD must be square, but it has 500 rows and 499 columns"),
    "LD's row and column names are the same" =
        list(d = local({
                 d <- D1
                 rownames(d$LD)[1] <- "x"
                 d
             }),
             error = "dataset 1: LD rownames != colnames"),
    "LD covers every snp" =
        list(d = local({
                 d <- D1
                 d$LD <- d$LD[-1, -1]
                 d
             }),
             error = "dataset 1: LD must cover every snp"),
    ## the values of the elements the first cases don't reach
    "pvalues are numeric" =
        list(d = local({
                 d <- pv1
                 d$pvalues <- as.character(d$pvalues)
                 d
             }),
             error = "dataset 1: pvalues must be numeric"),
    "pvalues are above 0 and at most 1" =
        list(d = local({
                 d <- pv1
                 d$pvalues[1] <- 0
                 d
             }),
             error = "dataset 1: pvalues must be above 0 and at most 1"),
    "beta is numeric" =
        list(d = local({
                 d <- D1
                 d$beta <- as.character(d$beta)
                 d
             }),
             error = "dataset 1: beta must be numeric"),
    "sdY is a single finite number above 0" =
        list(d = local({
                 d <- D1
                 d$sdY <- 0
                 d
             }),
             error = "dataset 1: sdY must be a single finite number above 0"),
    "position is above 0" =
        list(d = local({
                 d <- D1
                 d$position[1] <- 0L
                 d
             }),
             error = "dataset 1: position must hold finite whole numbers above 0"),
    "position holds whole numbers" =
        list(d = local({
                 d <- D1
                 d$position <- d$position + 0.5
                 d
             }),
             error = "dataset 1: position must hold finite whole numbers above 0"),
    "position is finite" =
        list(d = local({
                 d <- D1
                 d$position <- as.numeric(d$position)
                 d$position[1] <- Inf
                 d
             }),
             error = "dataset 1: position must hold finite whole numbers above 0"),
    "position is numeric" =
        list(d = local({
                 d <- D1
                 d$position <- as.character(d$position)
                 d
             }),
             error = "dataset 1: position must hold finite whole numbers above 0"),
    "position has no missing values" =
        list(d = local({
                 d <- D1
                 d$position[1] <- NA
                 d
             }),
             error = "dataset 1: position contains missing values"),
    "a whole-number position stored as a double passes" =
        list(d = local({
                 d <- D1
                 d$position <- as.numeric(d$position) * 1e6
                 d
             })),
    ## with no engine, every element given is checked, z apart
    "with no engine, pvalues are checked beside beta and varbeta" =
        list(d = local({
                 d <- D1
                 d$pvalues <- rep(2, length(d$snp))
                 d
             }),
             error = "dataset 1: pvalues must be above 0 and at most 1"),
    "with no engine, a case-control trait's sdY is checked" =
        list(d = local({
                 d <- D1
                 d$type <- "cc"
                 d$sdY <- 0
                 d
             }),
             error = "dataset 1: sdY must be a single finite number above 0"),
    "with no engine, z is not checked" =
        list(d = local({
                 d <- D1
                 d$z <- rep(NA_real_, length(d$snp))
                 d
             })),
    "with no engine, a dataset with only what SuSiE reads passes" =
        list(d = z1),
    "with no engine, the warning for a dataset without beta or p-values comes from z" =
        list(d = local({
                 d <- z1
                 d$z[] <- 0.1
                 d
             }),
             warning = "dataset 1: the smallest p-value is 0.92034, so this dataset may have no association in the region."),
    ## engine = "wakefield" checks what Wakefield reads
    "wakefield: a valid dataset passes" =
        list(args = wakefield, d = D1),
    "wakefield: a dataset given as p-values passes" =
        list(args = wakefield, d = pv1),
    "wakefield: a case-control trait with beta and varbeta needs neither s nor sdY" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$type <- "cc"
                 d$sdY <- NULL
                 d
             })),
    "wakefield: type is needed" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$type <- NULL
                 d
             }),
             error = "dataset 1: type is missing; set type to \"quant\" for a quantitative trait or \"cc\" for a case-control trait"),
    "wakefield: type is quant or cc" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$type <- "binary"
                 d
             }),
             error = "dataset 1: type must be quant or cc"),
    "wakefield: LD is not read" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$LD <- d$LD[, -1]
                 d
             })),
    "wakefield: pvalues beside beta and varbeta are not read" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$pvalues <- rep(2, length(d$snp))
                 d
             })),
    "wakefield: a case-control trait's sdY is not read" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$type <- "cc"
                 d$sdY <- 0
                 d
             })),
    "wakefield: MAF beside sdY is not read" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$MAF[1] <- 0
                 d
             })),
    "wakefield: without sdY, MAF is read" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$sdY <- NULL
                 d$MAF[1] <- 0
                 d
             }),
             error = "dataset 1: MAF should be a numeric, strictly >0 & <1"),
    "wakefield: beta is read" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$beta[1] <- NA
                 d
             }),
             error = "dataset 1: beta contains missing values"),
    "wakefield: the structure is still checked" =
        list(args = wakefield, d = local({
                 d <- D1
                 d$snp[2] <- d$snp[1]
                 d
             }),
             error = "dataset 1: duplicated snps found"),
    ## engine = "susie" checks what SuSiE reads
    "susie: a dataset with beta, varbeta, LD and N passes" =
        list(args = susie, d = D1),
    "susie: a dataset with z, LD and N passes" =
        list(args = susie, d = z1),
    "susie: LD is needed" =
        list(args = susie, d = local({
                 d <- D1
                 d$LD <- NULL
                 d
             }),
             error = "dataset 1: LD is missing; SuSiE needs LD, the correlation matrix of the snps"),
    "susie: N is needed" =
        list(args = susie, d = local({
                 d <- D1
                 d$N <- NULL
                 d
             }),
             error = "dataset 1: N is missing; SuSiE needs N, the sample size"),
    "susie: without z, beta and varbeta are needed" =
        list(args = susie, d = local({
                 d <- D1
                 d$beta <- NULL
                 d
             }),
             error = "dataset 1: beta is missing; without z, SuSiE needs beta and varbeta"),
    "susie: without z, beta or varbeta" =
        list(args = susie, d = local({
                 d <- z1
                 d$z <- NULL
                 d
             }),
             error = "dataset 1: beta and varbeta are missing; without z, SuSiE needs beta and varbeta"),
    "susie: z has no missing values" =
        list(args = susie, d = local({
                 d <- z1
                 d$z[1] <- NA
                 d
             }),
             error = "dataset 1: z must be a finite number for each snp"),
    "susie: z is finite" =
        list(args = susie, d = local({
                 d <- z1
                 d$z[1] <- Inf
                 d
             }),
             error = "dataset 1: z must be a finite number for each snp"),
    "susie: z is numeric" =
        list(args = susie, d = local({
                 d <- z1
                 d$z <- as.character(d$z)
                 d
             }),
             error = "dataset 1: z must be a finite number for each snp"),
    "susie: z has a value for each snp" =
        list(args = susie, d = local({
                 d <- z1
                 d$z <- d$z[-1]
                 d
             }),
             error = "dataset 1: z must be a finite number for each snp"),
    "susie: with z, beta is not read" =
        list(args = susie, d = local({
                 d <- D1
                 d$z <- d$beta / sqrt(d$varbeta)
                 d$beta[1] <- NA
                 d
             })),
    "susie: type is not read" =
        list(args = susie, d = local({
                 d <- D1
                 d$type <- "binary"
                 d
             })),
    "susie: LD is read" =
        list(args = susie, d = local({
                 d <- D1
                 d$LD <- d$LD * 1.5
                 d
             }),
             error = "dataset 1: LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal"),
    "susie: N is read" =
        list(args = susie, d = local({
                 d <- D1
                 d$N <- 1
                 d
             }),
             error = "dataset 1: N must be finite and > 1"),
    ## the arguments
    "warn.minp sets when the warning is given" =
        list(d = D1, args = list(warn.minp = 1e-30),
             warning = "dataset 1: the smallest p-value is < 2.22e-16, so this dataset may have no association in the region."),
    "req sets the required elements" =
        list(d = D1, args = list(req = c("snp", "pvalues")),
             error = "dataset 1: missing required element(s) pvalues"),
    ## datasets that fail two checks
    "a wrong type and duplicated snps" =
        list(d = local({
                 d <- D1
                 d$type <- "binary"
                 d$snp[2] <- d$snp[1]
                 d
             }),
             ## the structure is checked before any element's values
             error = "dataset 1: duplicated snps found"),
    "an s out of range, without beta or p-values" =
        list(d = local({
                 d <- D1
                 d$beta <- NULL
                 d$varbeta <- NULL
                 d$type <- "cc"
                 d$s <- 1.5
                 d
             }),
             error = "dataset 1: s must be between 0 and 1"),
    "a bad LD stops a weak dataset before the p-value warning" =
        list(d = local({
                 d <- D1
                 d$beta <- d$beta / 100
                 d$LD <- d$LD[, -1]
                 d
             }),
             error = "dataset 1: LD must be square, but it has 500 rows and 499 columns"))

for (what in names(cases)) {
    test_that(paste0("check_dataset(): ", what), {
        case <- cases[[what]]
        for (suffix in list(1, "")) {
            check <- function() do.call(check_dataset, c(list(case$d, suffix), case$args))
            ## without a label a message has no "dataset 1: ", and the
            ## not-a-list message begins "the dataset is"
            expected <- function(m) {
                if (identical(suffix, 1))
                    return(m)
                sub("^dataset 1: ", "", sub("^dataset 1: this is", "the dataset is", m))
            }
            ## a warning comes before any refusal
            until_refused <- function() tryCatch(check(), error = function(e) NULL)
            if (is.null(case$warning)) {
                expect_no_warning(until_refused())
            } else {
                w <- expect_warning(until_refused(), class = "coloc_min_p")
                expect_identical(sub("\n.*", "", conditionMessage(w)), expected(case$warning))
            }
            if (!is.null(case$error)) {
                e <- expect_error(suppressWarnings(check()))
                expect_identical(conditionMessage(e), expected(case$error))
            } else {
                expect_null(suppressWarnings(check()))
            }
        }
    })
}

## ---- what else check_dataset() must keep --------------------------------------
## check_dataset()'s errors and warnings show no call, as the message names the
## dataset

test_that("check_dataset()'s errors and warnings show no call", {
    dup <- D1
    dup$snp[2] <- dup$snp[1]
    weak <- D1
    weak$beta <- weak$beta / 100
    per_snp_N <- D1
    per_snp_N$N <- rep(D1$N, length(D1$snp))
    e <- expect_error(check_dataset(dup, 1), "dataset 1: duplicated snps found",
                      fixed = TRUE)
    expect_null(conditionCall(e))
    w <- expect_warning(check_dataset(weak, 1),
                        "dataset 1: the smallest p-value is 0.93383, so this dataset may have no association in the region.",
                        fixed = TRUE, class = "coloc_min_p")
    expect_null(conditionCall(w))
    e <- expect_error(check_dataset(1:3, 1),
                      "dataset 1: this is a numeric vector, not a list of coloc_dataset elements (snp, beta, varbeta, ...)",
                      fixed = TRUE)
    expect_null(conditionCall(e))
    e <- expect_error(check_dataset(per_snp_N, 1),
                      "dataset 1: N must be a single number; coloc does not take per-snp sample sizes",
                      fixed = TRUE)
    expect_null(conditionCall(e))
})

test_that("check_dataset() refuses a data.frame, as coloc() does", {
    df <- data.frame(snp = D1$snp, beta = D1$beta, varbeta = D1$varbeta,
                     type = "quant", sdY = D1$sdY)
    e <- expect_error(check_dataset(df))
    expect_identical(conditionMessage(e),
                     "the dataset is a data.frame, not a list of coloc_dataset elements (snp, beta, varbeta, ...)")
    ## one row, so type is a single value and only the data.frame test refuses it
    e <- expect_error(check_dataset(df[1, ], 1))
    expect_identical(conditionMessage(e),
                     "dataset 1: this is a data.frame, not a list of coloc_dataset elements (snp, beta, varbeta, ...)")
})

test_that("engine must be NULL, wakefield or susie", {
    for (engine in list("abf", c("wakefield", "susie"), 1, NA))
        expect_error(check_dataset(D1, engine = engine),
                     "engine must be NULL or one of: wakefield, susie", fixed = TRUE)
})

test_that("a check for Wakefield doesn't test LD's values, and one for SuSiE does", {
    ## an LD that is no correlation matrix, in a dataset already made
    D <- as_coloc_dataset(D1)
    D$LD <- D$LD * 1.5
    expect_null(check_dataset(D, engine = "wakefield"))
    expect_error(check_dataset(D, engine = "susie"), "LD must be a correlation matrix", fixed = TRUE)
    expect_error(check_dataset(D), "LD must be a correlation matrix", fixed = TRUE)
})
