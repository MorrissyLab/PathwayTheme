#!/usr/bin/env Rscript
# Reference limma run, used to validate PathwayTheme's pure-Python port of
# limma's empirical-Bayes moderated t (fit_fdist / trigamma_inverse).
#
# Reads a pathway score matrix and a two-group design, fits exactly the model
# the Python side fits -- ~ 0 + group, contrast case - reference -- and writes
# both the per-pathway statistics AND the fitted hyperparameters (d0, s0^2),
# so the port can be checked at the prior level, not just at the t-statistic
# level where errors could cancel.
#
#   Rscript limma_reference.R <scores.tsv> <groups.tsv> <case> <ref> <out.tsv>
#
# scores.tsv : pathways x samples, tab separated, row names in column 1
# groups.tsv : two columns, "sample" and "group"

suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("need 5 arguments")
scores_path <- args[1]; groups_path <- args[2]
case <- args[3]; ref <- args[4]; out_path <- args[5]

scores <- as.matrix(read.delim(scores_path, row.names = 1, check.names = FALSE))
groups <- read.delim(groups_path, stringsAsFactors = FALSE)
rownames(groups) <- groups$sample

keep <- groups$sample[groups$group %in% c(case, ref)]
keep <- intersect(keep, colnames(scores))
X <- scores[, keep, drop = FALSE]
g <- factor(groups[keep, "group"], levels = c(ref, case))

design <- model.matrix(~ 0 + g)
colnames(design) <- c(ref, case)

fit <- lmFit(X, design)
contrast <- makeContrasts(contrasts = paste0(case, " - ", ref),
                          levels = design)
fit2 <- contrasts.fit(fit, contrast)
eb <- eBayes(fit2)

res <- data.frame(
  pathway     = rownames(X),
  effect      = as.numeric(eb$coefficients[, 1]),
  s2_ordinary = as.numeric(eb$sigma^2),
  s2_post     = as.numeric(eb$s2.post),
  t_moderated = as.numeric(eb$t[, 1]),
  p_moderated = as.numeric(eb$p.value[, 1]),
  df_residual = as.numeric(eb$df.residual),
  df_prior    = as.numeric(eb$df.prior),
  s2_prior    = as.numeric(eb$s2.prior),
  stringsAsFactors = FALSE
)
write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("limma %s | %s vs %s | n=%d+%d | d0=%.6f s0^2=%.6f\n",
            as.character(packageVersion("limma")), case, ref,
            sum(g == case), sum(g == ref),
            eb$df.prior[1], eb$s2.prior[1]))
