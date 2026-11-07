# Precompile the validation vignette
# =================================
#
# The vignette compares against the New York Fed's published ACM series,
# which the package downloads but does not redistribute (see LICENSE.note).
# A vignette that downloads cannot be rebuilt on CRAN, so the source,
# validation.Rmd.orig, is knitted here into validation.Rmd: the code and its
# output as static text, plus the figures in figures/. R CMD build then only
# renders that markdown, and needs neither the data nor a network.
#
# Edit validation.Rmd.orig, never validation.Rmd, then re-run this script.
# Re-run it too after any change that could move the vignette's numbers, and
# read the prose against the new output: several of them are quoted in the
# text.
#
# Needs an internet connection and the readxl package.
#
# Usage, from the package root: Rscript vignettes/precompile.R
# =========================================================================

stopifnot(file.exists("DESCRIPTION"),
          requireNamespace("knitr", quietly = TRUE),
          requireNamespace("readxl", quietly = TRUE))

# Knit against this source tree, not whatever version happens to be installed.
lib <- tempfile("termpremia-lib-")
dir.create(lib)
install.packages(".", repos = NULL, type = "source", lib = lib, quiet = TRUE)
.libPaths(c(lib, .libPaths()))

old <- setwd("vignettes")
unlink("figures", recursive = TRUE)
knitr::knit("validation.Rmd.orig", output = "validation.Rmd", quiet = TRUE)
setwd(old)

message("Wrote vignettes/validation.Rmd and vignettes/figures/.")
