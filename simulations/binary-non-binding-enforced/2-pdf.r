# One PDF per method (pdf/initial-<method>.pdf) with a summary page and a
# sample-size plot per config.
# Usage: Rscript 2-pdf.r 1 in out pdf

args <- commandArgs(trailingOnly = TRUE)

library(MAMS)
library(MAMSS)
library(dotfunctions)
options(dotfunctions_dontask = TRUE)

path_in  <- paste0(args[2], "/")
path_out <- paste0(args[3], "/")
path_pdf <- paste0(args[4], "/")

print(load(paste0(path_in, "par")))

# power, FWER and expected sample sizes of one MAMSS result
metrics <- function(m) {
  pow  <- m$H1$target$global$efficacy[1]   # P(reject >= 1 arm | H1) = power
  fwer <- m$H0$target$global$efficacy[1]   # P(reject >= 1 arm | H0) = FWER

  # drop the last two book-keeping columns (see plot.mamss)
  s1   <- m$H1$sample[, -(ncol(m$H1$sample) + (-1:0)), drop = FALSE]
  s0   <- m$H0$sample[, -(ncol(m$H0$sample) + (-1:0)), drop = FALSE]
  list(
    power   = pow,
    fwer    = fwer,
    ess.h1  = mean(rowSums(s1)),   # expected sample size under H1
    ess.h0  = mean(rowSums(s0)),   # expected sample size under H0
    maxn    = max(m$look$n)        # planned maximum sample size
  )
}

# plain-text PDF page
text_page <- function(lines, cex = 0.8) {
  plot.new()
  par(mar = c(0, 0, 0, 0), family = "mono")
  n <- length(lines)
  y <- seq(0.97, 0.03, length.out = max(n, 1))
  text(0.02, y, lines, adj = c(0, 0.5), cex = cex)
  par(family = "")
}

dir.create(path_pdf, showWarnings = FALSE, recursive = TRUE)

ids     <- seq(.an(args[1]), n.config, 1)
methods <- unique(id.config$method[ids])

for (meth in methods) {
  pdf_file <- paste0(path_pdf, "initial-", meth, ".pdf")
  pdf(pdf_file, width = 8.5, height = 11, onefile = TRUE)

  for (cw in ids[id.config$method[ids] == meth]) {
    pos  <- id.config$pos[cw]
    file <- paste0(path_out, "initial/", pos)
    if (!file.exists(file)) {
      cat("skip config", pos, "- no result file\n")
      next
    }
    load(file)                       # provides 'out'
    m  <- out$mamss0.ptest
    mm <- metrics(m)
    pc <- out$rates["pc"]
    pt <- out$rates["pt"]

    cat("page for config", pos, "(", meth, ")\n")

    header <- c(
      sprintf("Config %s   (%s)   [BINARY, NON-BINDING, FUT ENFORCED]", pos, id.config$id[cw]),
      strrep("-", 60),
      sprintf("  K (treatment arms)  : %s", id.config$K[cw]),
      sprintf("  J (looks)           : %s", id.config$J[cw]),
      sprintf("  p  (target win-prob) : %s", id.config$p[cw]),
      sprintf("  p0 (null win-prob)   : %s", id.config$p0[cw]),
      sprintf("  method               : %s", id.config$method[cw]),
      "",
      sprintf("  true response rates  : control pi_c = %.3f", pc),
      sprintf("                         effect  pi_t = %.3f  (= pc + 2p-1)", pt),
      "",
      sprintf("  target alpha = %-5s   target power = %s",
              id.config$alpha[cw], id.config$power[cw]),
      strrep("-", 60),
      "  SIMULATED OPERATING CHARACTERISTICS  (Bernoulli, stats.ptest)",
      sprintf("    Power  (P reject >=1 | H1) : %.4f", mm$power),
      sprintf("    FWER   (P reject >=1 | H0) : %.4f", mm$fwer),
      sprintf("    Max sample size (planned)  : %.0f", mm$maxn),
      sprintf("    ESS under H1               : %.1f", mm$ess.h1),
      sprintf("    ESS under H0               : %.1f", mm$ess.h0),
      "",
      "  NB: non-binding (wider) boundary with futility ENFORCED, so power",
      "      and ESS are realistic. Design powered for the binary test",
      "      (delta/sd); residual gap from target is small-sample discreteness."
    )
    text_page(header)

    plot(m, type = "size", hypothesis = "H1",
         title = sprintf("Config %s - sample size (H1) [binary]", pos))
  }

  dev.off()
  cat("\n\t Wrote", pdf_file, "\n")
}
if (!interactive()) q("no")
