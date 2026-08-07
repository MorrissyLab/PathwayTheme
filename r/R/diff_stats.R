## Per-pathway two-group statistics used by the differential stage.
##
## Three methods, all operating on one pathway x sample score matrix split into
## a *case* and a *reference* set of columns:
##
## * `welch`        -- Welch's unequal-variance t-test.
## * `mannwhitney`  -- Mann-Whitney U rank test.
## * `moderated_t`  -- limma's empirical-Bayes moderated t-statistic
##   (Smyth 2004).  [pt_fit_fdist()] and [pt_trigamma_inverse()] mirror
##   `limma::fitFDist` / `limma::trigammaInverse`; the two-group case reduces to
##   a pooled linear model with variance shrinkage.
##
## Each method returns three equal-length vectors over the pathways (rows):
## `effect` (mean_case - mean_reference), `statistic`, `p_value`.

# ---- limma variance-model helpers -----------------------------------------

#' Solve `trigamma(y) = x` for `y`
#'
#' A port of `limma::trigammaInverse`.
#'
#' @param x a numeric scalar.
#' @return The solution `y`.
#' @export
pt_trigamma_inverse <- function(x) {
  x <- as.double(x)
  if (!is.finite(x)) return(NA_real_)
  if (x > 1e7) return(1 / sqrt(x))
  if (x < 1e-6) return(1 / x)
  y <- 0.5 + 1 / x
  for (i in seq_len(50L)) {
    tri <- trigamma(y)
    dif <- tri * (1 - tri / x) / psigamma(y, deriv = 2L)
    y <- y + dif
    if (-dif / y < 1e-8) break
  }
  y
}

#' Estimate the scaled-F prior for a set of sample variances
#'
#' A port of `limma::fitFDist`.  Given per-pathway residual variances `var`,
#' each on `df1` degrees of freedom, returns the prior degrees of freedom and
#' prior variance used for empirical-Bayes shrinkage.
#'
#' @param var numeric vector of residual variances.
#' @param df1 numeric vector (or scalar) of residual degrees of freedom.
#' @return A list with `d0` and `s0_sq`.
#' @export
pt_fit_fdist <- function(var, df1) {
  var <- as.double(var)
  df1 <- rep_len(as.double(df1), length(var))

  ok <- is.finite(var) & var > 0 & is.finite(df1) & df1 > 0
  x <- var[ok]
  g1 <- df1[ok]
  if (!length(x)) return(list(d0 = Inf, s0_sq = NA_real_))

  z <- log(x)
  e <- z - digamma(g1 / 2) + log(g1 / 2)
  emean <- mean(e)
  n <- length(e)
  if (n <= 1L) return(list(d0 = Inf, s0_sq = exp(emean)))

  evar <- sum((e - emean)^2) / (n - 1)
  evar <- evar - mean(trigamma(g1 / 2))

  if (evar > 0) {
    d0 <- 2 * pt_trigamma_inverse(evar)
    s0_sq <- exp(emean + digamma(d0 / 2) - log(d0 / 2))
  } else {
    d0 <- Inf
    s0_sq <- exp(emean)
  }
  list(d0 = d0, s0_sq = s0_sq)
}


# ---- the three tests -------------------------------------------------------

#' Two-group tests over a pathway x sample matrix
#'
#' @param case,ref numeric matrices with the same rows (pathways) and the
#'   samples of each group in the columns.
#' @return A list of three vectors: `effect`, `statistic`, `p_value`.
#' @export
pt_welch <- function(case, ref) {
  n1 <- rowSums(!is.na(case))
  n2 <- rowSums(!is.na(ref))
  m1 <- rowMeans(case, na.rm = TRUE)
  m2 <- rowMeans(ref, na.rm = TRUE)
  v1 <- .pt_row_var(case, m1, n1)
  v2 <- .pt_row_var(ref, m2, n2)

  se1 <- v1 / n1
  se2 <- v2 / n2
  denom <- sqrt(se1 + se2)
  stat <- (m1 - m2) / denom
  df <- (se1 + se2)^2 / (se1^2 / (n1 - 1) + se2^2 / (n2 - 1))
  pval <- rep(NA_real_, length(stat))
  fin <- is.finite(stat) & is.finite(df) & df > 0
  pval[fin] <- 2 * stats::pt(-abs(stat[fin]), df = df[fin])
  list(effect = m1 - m2, statistic = stat, p_value = pval)
}

#' @rdname pt_welch
#' @export
pt_mannwhitney <- function(case, ref) {
  n <- nrow(case)
  effect <- rowMeans(case, na.rm = TRUE) - rowMeans(ref, na.rm = TRUE)
  stat <- rep(NA_real_, n)
  pval <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    x <- case[i, ][!is.na(case[i, ])]
    y <- ref[i, ][!is.na(ref[i, ])]
    if (!length(x) || !length(y)) next
    # scipy's method='auto': exact when there are no ties and the smaller group
    # has at most 8 observations, asymptotic (tie- and continuity-corrected)
    # otherwise.  R's wilcox.test implements both.
    use_exact <- !anyDuplicated(c(x, y)) && min(length(x), length(y)) <= 8L
    res <- tryCatch(
      suppressWarnings(stats::wilcox.test(x, y, alternative = "two.sided",
                                          exact = use_exact, correct = TRUE)),
      error = function(e) NULL)
    if (is.null(res)) next
    stat[[i]] <- as.double(res$statistic)
    pval[[i]] <- as.double(res$p.value)
  }
  list(effect = effect, statistic = stat, p_value = pval)
}

#' @rdname pt_welch
#' @export
pt_moderated_t <- function(case, ref) {
  n1 <- rowSums(!is.na(case))
  n2 <- rowSums(!is.na(ref))
  m1 <- rowMeans(case, na.rm = TRUE)
  m2 <- rowMeans(ref, na.rm = TRUE)
  effect <- m1 - m2

  ss1 <- rowSums((case - m1)^2, na.rm = TRUE)
  ss2 <- rowSums((ref - m2)^2, na.rm = TRUE)
  dfr <- n1 + n2 - 2                              # residual df per pathway
  s2 <- ss1 + ss2
  s2 <- ifelse(dfr > 0, s2 / dfr, NA_real_)
  stdev_unscaled <- sqrt(1 / n1 + 1 / n2)

  good <- is.finite(s2) & dfr > 0
  prior <- pt_fit_fdist(s2[good], dfr[good])
  d0 <- prior$d0
  s0_sq <- prior$s0_sq

  if (!is.finite(d0)) {
    s2_post <- rep(s0_sq, length(s2))
    df_total <- rep(Inf, length(dfr))
  } else {
    s2_post <- (d0 * s0_sq + dfr * s2) / (d0 + dfr)
    df_total <- dfr + d0
  }

  stat <- effect / (stdev_unscaled * sqrt(s2_post))
  pval <- rep(NA_real_, length(stat))
  fin <- is.finite(stat)
  inf_df <- fin & !is.finite(df_total)
  fin_df <- fin & is.finite(df_total) & df_total > 0
  pval[fin_df] <- 2 * stats::pt(-abs(stat[fin_df]), df = df_total[fin_df])
  pval[inf_df] <- 2 * stats::pnorm(-abs(stat[inf_df]))
  list(effect = effect, statistic = stat, p_value = pval)
}

#' Dispatch a two-group test by name
#'
#' @param method one of `"welch"`, `"mannwhitney"`, `"moderated_t"`.
#' @param case,ref pathway x sample matrices for the two groups.
#' @return A list of `effect`, `statistic`, `p_value`.
#' @export
pt_run_test <- function(method, case, ref) {
  fn <- switch(method,
               welch = pt_welch,
               mannwhitney = pt_mannwhitney,
               moderated_t = pt_moderated_t,
               stop("Unknown diff method '", method,
                    "'; expected one of mannwhitney, moderated_t, welch",
                    call. = FALSE))
  fn(case, ref)
}

#' Benjamini-Hochberg adjusted p-values
#'
#' `NA`s pass through and are excluded from the ranking, matching the Python
#' implementation.
#'
#' @param p numeric vector of p-values.
#' @return A numeric vector of the same length.
#' @export
pt_benjamini_hochberg <- function(p) {
  p <- as.double(p)
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  m <- sum(ok)
  if (m == 0L) return(out)
  idx <- which(ok)
  ord <- idx[order(p[idx], method = "radix")]
  ranked <- p[ord] * m / seq_len(m)
  ranked <- rev(cummin(rev(ranked)))          # enforce monotonicity
  out[ord] <- pmin(pmax(ranked, 0), 1)
  out
}

## Row-wise sample variance (ddof = 1) given precomputed means and counts.
.pt_row_var <- function(m, means, n) {
  ss <- rowSums((m - means)^2, na.rm = TRUE)
  ifelse(n > 1L, ss / (n - 1), NA_real_)
}
