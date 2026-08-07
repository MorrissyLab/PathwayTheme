#!/usr/bin/env Rscript
## Launcher for the pathwaytheme command-line interface.
suppressPackageStartupMessages(library(pathwaytheme))
status <- pt_main()
quit(status = if (is.numeric(status)) as.integer(status) else 0L)
