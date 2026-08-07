## Command-line interface: `pathwaytheme run config.yaml`.
##
## Call it directly, or through the launcher installed at
## `system.file("scripts", "pathwaytheme.R", package = "pathwaytheme")`:
##
##     Rscript -e 'pathwaytheme::pt_main()' run config.yaml
##     Rscript "$(Rscript -e 'cat(system.file("scripts/pathwaytheme.R",
##                                            package="pathwaytheme"))')" \
##             run config.yaml

#' Command-line entry point
#'
#' Two subcommands, matching the Python CLI:
#'
#' * `run <config.yaml> [--output-dir DIR] [--no-figures] [--quiet]`
#' * `init [-o FILE]` -- print or write a default configuration
#'
#' @param argv character vector of arguments; defaults to the process's
#'   trailing command-line arguments.
#' @return `0` on success, invisibly; other values signal an error.
#' @export
pt_main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (!length(argv) || argv[[1L]] %in% c("-h", "--help", "help"))
    return(.pt_cli_usage())
  command <- argv[[1L]]
  rest <- argv[-1L]
  if (identical(command, "run")) return(.pt_cli_run(rest))
  if (identical(command, "init")) return(.pt_cli_init(rest))
  message("unknown command '", command, "'")
  .pt_cli_usage()
  invisible(2L)
}

.pt_cli_usage <- function() {
  cat(paste(
    "pathwaytheme -- omic enrichment (ssGSEA/EnrichR/GoSlim) + PCA pathway analysis",
    "",
    "usage:",
    "  pathwaytheme run <config.yaml> [--output-dir DIR] [--no-figures] [--quiet]",
    "  pathwaytheme init [-o FILE]",
    "", sep = "\n"))
  invisible(0L)
}

.pt_cli_flag <- function(args, flag) flag %in% args

.pt_cli_value <- function(args, flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1L]]
}

.pt_cli_run <- function(args) {
  positional <- args[!startsWith(args, "-")]
  # drop the value that follows --output-dir
  od <- .pt_cli_value(args, "--output-dir")
  if (!is.null(od)) positional <- setdiff(positional, od)
  if (!length(positional)) {
    message("run: a config YAML path is required")
    return(invisible(2L))
  }
  config <- pt_config_from_yaml(positional[[1L]])
  if (!is.null(od)) config$output_dir <- od
  if (.pt_cli_flag(args, "--no-figures")) config$viz$make_figures <- FALSE
  pt_run(config, verbose = !.pt_cli_flag(args, "--quiet"))
  invisible(0L)
}

.pt_cli_init <- function(args) {
  cfg <- pt_config()
  out <- .pt_cli_value(args, "-o", .pt_cli_value(args, "--output"))
  if (!is.null(out)) {
    pt_config_to_yaml(cfg, out)
    cat("wrote default config to ", out, "\n", sep = "")
  } else {
    cat(yaml::as.yaml(unclass(cfg), handlers = list(
      `NULL` = function(x) structure("~", class = "verbatim"))))
  }
  invisible(0L)
}
