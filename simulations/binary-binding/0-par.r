# Parameter grid for binary-binding: writes in/par and resets out/.
# Usage: Rscript 0-par.r
#
# p and p0 are MAMS rank effects. With control response rate pc the true
# rates are pi_c = pc and pi_t = pc + (2p - 1); only the (pc, p) pairs in
# valid_pairs are run.

rm(list = ls())

# dotfunctions: github.com/dlc48/dotfunctions
if (!requireNamespace("dotfunctions", quietly = TRUE)) {
  devtools::install_github("dlc48/dotfunctions")
}
library(dotfunctions)
options(dotfunctions_dontask = TRUE)

.idf(c(1, 4), "K")
.idf(c(1, 2), "J")
.idf(c(0.52, 0.55, 0.60, 0.65, 0.75, 0.85, 0.95), "p")
.idf(c(0.5), "p0")
.idf(c(0.05, 0.10, 0.30, 0.50, 0.70, 0.90), "pc")
.idf(c(0.05), "alpha")
.idf(c(0.9), "power")
.idf(c("simultaneous", "sep"), "method")

# all combinations of the parameter sets
.idc("config", id.alpha, id.power, id.K, id.J, id.p, id.p0, id.pc, id.method)

id.config$K      <- .anac(id.config$K)
id.config$J      <- .anac(id.config$J)
id.config$p      <- .anac(id.config$p)
id.config$p0     <- .anac(id.config$p0)
id.config$pc     <- .anac(id.config$pc)
id.config$alpha  <- .anac(id.config$alpha)
id.config$power  <- .anac(id.config$power)
id.config$method <- as.character(id.config$method)

valid_pairs <- data.frame(
  pc = c(0.05, 0.05, 0.10, 0.30, 0.30, 0.50, 0.50, 0.70, 0.70, 0.90),
  p  = c(0.65, 0.95, 0.85, 0.65, 0.75, 0.65, 0.55, 0.60, 0.52, 0.52)
)

id.config <- id.config[paste(id.config$pc, id.config$p) %in%
                          paste(valid_pairs$pc, valid_pairs$p), ]
n.config  <- nrow(id.config)

# reset output folder
system("rm -rf out/*")
system("mkdir -p out/initial")

# save grid for 1-run.r
save(list = c(.ide()), file = "in/par")
