# Parameter grid for non-binding-enforced: writes in/par and resets out/.
# Usage: Rscript 0-par.r

rm(list = ls())

# dotfunctions: github.com/dlc48/dotfunctions
if (!requireNamespace("dotfunctions", quietly = TRUE)) {
  devtools::install_github("dlc48/dotfunctions")
}
library(dotfunctions)
options(dotfunctions_dontask = TRUE)

.idf(c(1, 4), "K")
.idf(c(1, 2), "J")
.idf(c(0.65, 0.75), "p")
.idf(c(0.5), "p0")
.idf(c(0.05), "alpha")
.idf(c(0.9), "power")
.idf(c("simultaneous", "sep"), "method")

# all combinations of the parameter sets
.idc("config", id.alpha, id.power, id.K, id.J, id.p, id.p0, id.method)

id.config$K      <- .anac(id.config$K)
id.config$J      <- .anac(id.config$J)
id.config$p      <- .anac(id.config$p)
id.config$p0     <- .anac(id.config$p0)
id.config$alpha  <- .anac(id.config$alpha)
id.config$power  <- .anac(id.config$power)
id.config$method <- as.character(id.config$method)

# reset output folder
system("rm -rf out/*")
system("mkdir -p out/initial")

# save grid for 1-run.r
save(list = c(.ide()), file = "in/par")
