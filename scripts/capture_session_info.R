#!/usr/bin/env Rscript

output_file <- if (length(commandArgs(trailingOnly = TRUE)) >= 1L) {
  commandArgs(trailingOnly = TRUE)[1L]
} else {
  file.path("results", "sessionInfo.txt")
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), output_file)
message("Saved session information to ", normalizePath(output_file))
