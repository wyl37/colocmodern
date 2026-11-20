##' Wakefield's approximate Bayes factor for each snp
##'
##' `wakefield_bf()` gives Wakefield's approximate Bayes factor (ABF) for each
##' snp in a dataset, on the log scale. Each is the evidence for an effect at
##' that snp against no effect. [finemap()] computes these with
##' `engine = "wakefield"`, and [coloc()] compares the two traits snp by snp.
##'
##' Wakefield showed we could get ABF as
##'
##' log ABF = (log(1 - r) + r z^2) / 2
##'
##' We assume the true effect size ~ N(0,W), and set r = W / (W + V)
##' where V is the variance of the estimated effect size. z is the
##' usual Z score (beta/se(beta)). W is given by the `effect_priors`
##' parameter.
##'
##' Where the dataset includes beta and varbeta, this is
##' straightforward. Where it does not, we approximate the result from
##' pvalues, MAF and N. Please give beta and varbeta where available.
##'
##' For a quantitative trait the effect prior is in units of the
##' trait's standard deviation, so it is multiplied by `sdY`; when the
##' dataset has no `sdY`, it is estimated from `varbeta`, `MAF` and
##' `N`, with a warning. Again, please give sdY if known.
##'
##' @param d a dataset made with [coloc_dataset()] or [as_coloc_dataset()]
##' @param effect_priors the prior standard deviation of a causal snp's
##'     effect, named by the type of trait: `quant` in units of the trait's
##'     standard deviation, and `cc` on the log odds ratio scale. This option
##'     exists for special cases - do not change the standard priors unless
##'     you know what you are doing!
##' @return a data frame with one row per snp and the columns `snp`,
##'     `position` (`NA` when the dataset has none), `z`, `V` (the variance
##'     of the effect estimate), `r` (the shrinkage factor) and `lABF` (the
##'     log ABF)
##' @references Wakefield J (2009). Bayes factors for genome-wide association
##'     studies: comparison with P-values. Genetic Epidemiology 33(1),
##'     79-86. \doi{10.1002/gepi.20359}
##' @seealso [finemap()], which puts these in a `coloc_finemap`;
##' [check_dataset()], which gives more detail on how to structure a coloc dataset
##' @author Wei-Yu Lin, Claudia Giambartolomei, Vincent Plagnol, Chris Wallace
##' @importFrom stats qnorm lm coef
##' @noRd
wakefield_bf <- function(d, effect_priors = c(quant = 0.15, cc = 0.2)) {
    .need_coloc_dataset(d, "d")
    if (!is.numeric(effect_priors) ||
        !all(c("quant", "cc") %in% names(effect_priors)) ||
        !all(is.finite(effect_priors[c("quant", "cc")]) &
             effect_priors[c("quant", "cc")] > 0))
        stop("effect_priors should be positive numbers named quant and cc",
             call. = FALSE)
    ## checked for what Wakefield reads: a coloc_dataset may lack it, and its
    ## elements may have been changed since it was made
    check_dataset(d, engine = "wakefield")
    nd <- names(d)
    ## by the type's label: [[ would take a factor type by its code
    sd_prior <- effect_priors[[if (d$type == "quant") "quant" else "cc"]]
    if ("beta" %in% nd && "varbeta" %in% nd) {
        z <- d$beta / sqrt(d$varbeta)
        V <- d$varbeta
        ## a quantitative trait's prior is in units of its sd, estimated, with
        ## a warning, when not given
        if (d$type == "quant")
            sd_prior <- sd_prior *
                (if ("sdY" %in% nd) d$sdY else sdY.est(d$varbeta, d$MAF, d$N))
    } else {
        z <- qnorm(0.5 * d$pvalues, lower.tail = FALSE)
        V <- if (d$type == "quant") Var.data(d$MAF, d$N)
             else Var.data.cc(d$MAF, d$N, d$s)
    }
    r <- sd_prior^2 / (sd_prior^2 + V)
    data.frame(snp = d$snp,
               position = if ("position" %in% nd) d$position
                          else rep(NA_integer_, length(z)),
               z = z, V = V, r = r, lABF = 0.5 * (log(1 - r) + (r * z^2)),
               row.names = NULL)
}

##' variance of MLE of beta for quantitative trait, assuming var(y)=1
##'
##' Internal function
##' @title Var.data
##' @param f minor allele freq
##' @param N sample number
##' @return variance of MLE beta
##' @author Claudia Giambartolomei
##' @noRd
Var.data <- function(f, N) {
    1 / (2 * N * f * (1 - f))
}

##' variance of MLE of beta for case-control
##'
##' Internal function
##' @title Var.data
##' @inheritParams Var.data
##' @param s proportion of samples that are cases
##' @return variance of MLE beta
##' @author Claudia Giambartolomei
##' @noRd
Var.data.cc <- function(f, N, s) {
    1 / (2 * N * f * (1 - f) * s * (1 - s))
}

##' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
##'
##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
##' var(X) = 2*maf*(1-maf)
##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
##'
##' @title Estimate trait variance, internal function
##' @param vbeta vector of variance of coefficients
##' @param maf vector of MAF (same length as vbeta)
##' @param n sample size
##' @return estimated standard deviation of Y
##'
##' @author Chris Wallace
##' @noRd
sdY.est <- function(vbeta, maf, n) {
    warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
    oneover <- 1/vbeta
    nvx <- 2 * n * maf * (1-maf)
    m <- lm(nvx ~ oneover - 1)
    cf <- coef(m)[['oneover']]
    if(cf < 0)
        stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
    return(sqrt(cf))
}
