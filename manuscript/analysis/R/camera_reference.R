#!/usr/bin/env Rscript
# limma::camera on the gene-level matrix, as the competitive test that carries
# the inter-gene correlation term PathwayTheme's score-level test does not.
#
# PathwayTheme scores each gene set per sample (ssGSEA) and then runs a
# moderated t on the score matrix, one row per set.  That treats the sets as
# independent features.  camera asks the competitive question directly on the
# genes and inflates the variance of each set's statistic by
# 1 + (m - 1) * inter.gene.cor, so the two differ in exactly the term at issue.
#
#   Rscript camera_reference.R <expr.tsv> <groups.tsv> <sets.gmt> <case> <ref> <out.tsv>
#
# expr.tsv   : genes x samples, tab separated, row names in column 1
# groups.tsv : two columns, "sample" and "group"
# sets.gmt   : GMT; only sets with >= MIN and <= MAX measured members are tested

suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) stop("need 6 arguments")
expr_path <- args[1]; groups_path <- args[2]; gmt_path <- args[3]
case <- args[4]; ref <- args[5]; out_path <- args[6]

MIN <- 10; MAX <- 500          # the size window the ssGSEA run used

expr <- as.matrix(read.delim(expr_path, row.names = 1, check.names = FALSE))
groups <- read.delim(groups_path, stringsAsFactors = FALSE)
rownames(groups) <- groups$sample

keep <- intersect(groups$sample[groups$group %in% c(case, ref)], colnames(expr))
y <- expr[, keep, drop = FALSE]
g <- factor(groups[keep, "group"], levels = c(ref, case))

design <- model.matrix(~ 0 + g)
colnames(design) <- c(ref, case)
contrast <- makeContrasts(contrasts = paste0(case, " - ", ref), levels = design)

# ---- gene sets, as row indices into y --------------------------------------
lines <- readLines(gmt_path, warn = FALSE)
sets <- lapply(strsplit(lines, "\t"), function(p) p[-c(1, 2)])
names(sets) <- vapply(strsplit(lines, "\t"), `[`, "", 1)
index <- ids2indices(sets, rownames(y), remove.empty = TRUE)
sizes <- lengths(index)
index <- index[sizes >= MIN & sizes <= MAX]

res <- camera(y, index, design, contrast = contrast, inter.gene.cor = 0.01)
res <- data.frame(pathway = rownames(res), res, row.names = NULL,
                  stringsAsFactors = FALSE)
write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("camera %s | %s vs %s | n=%d+%d | %d sets | %d at FDR<0.05\n",
            as.character(packageVersion("limma")), case, ref,
            sum(g == case), sum(g == ref), nrow(res),
            sum(res$FDR < 0.05)))
