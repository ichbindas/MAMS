# Simulate every config of the grid: MAMS design, MAMSS evaluation.
# Usage: Rscript 1-run.r 1 in out
#
# Binding futility: arms crossing the lower boundary are dropped.

args <- commandArgs(trailingOnly = TRUE)

library(MAMS)
library(MAMSS)
library(dotfunctions)
options(dotfunctions_dontask = TRUE)

# load the grid written by 0-par.r
path_in  <- paste0(args[2], "/")
path_out <- paste0(args[3], "/")

print(load(paste0(path_in, "par")))
cat("\n\t", date(), "\n")

# process configs from args[1] onwards
id.configw <- id.config[seq(.an(args[1]), n.config, 1), ]
n.configw  <- nrow(id.configw)

for (cw in 1:n.configw) {
  cat("START config", id.configw$id[cw], "(", cw, "/", n.configw, ")\n")

  # skip configs that already have a result
  if (!any(dir(paste0(path_out, "initial/")) == id.configw$pos[cw])) {
    K <- id.configw$K[cw]
    J <- id.configw$J[cw]

    cat("\tDesign (MAMS::mams)\n")

    # 1. design
    mams0 <- MAMS::mams(
      K          = K,
      J          = J,
      alpha      = id.configw$alpha[cw],
      power      = id.configw$power[cw],
      p          = id.configw$p[cw],
      p0         = id.configw$p0[cw],
      r          = 1:J,
      r0         = 1:J,
      ushape     = "obf",
      lshape     = "triangular",
      method     = id.configw$method[cw],
      binding    = TRUE,
      sample.size = TRUE,
      nsim       = NULL,
      parallel   = FALSE,
      print      = TRUE
    )

    # 2. per-stage (incremental) sample sizes for MAMSS
    if (J > 1) {
      mMat <- cbind(
        mams0$rMat[, 1],
        mams0$rMat[, 2:J] - mams0$rMat[, 1:(J - 1)]
      ) * mams0$n
    } else {
      mMat <- mams0$rMat * mams0$n
    }

    # 3. operating characteristics with MAMSS: T1 has effect
    #    sqrt(2) * qnorm(p), the other arms none; sd = exp(0) = 1 (log link)
    delta  <- sqrt(2) * qnorm(id.configw$p[cw])
    coef.mu <- c(0, delta, rep(0, K - 1))

    coef.sigma <- 0

    # sep: arms stop individually, the trial never stops as a whole
    if (id.configw$method[cw] == "sep") {
      eff_trial <- function(eff.target) FALSE
      fut_trial <- function(fut.target) FALSE
    } else {
      eff_trial <- eff.trial.any
      fut_trial <- fut.trial.all
    }

    cat("\tEvaluation (MAMSS::mamss - stats.ztest)\n")

    mamss0.ztest <- mamss(
      model           = outcome ~ intervention,
      var             = list(outcome   = rNO,
                             intervention = alloc.balanced),
      coef.mu         = coef.mu,
      coef.sigma      = coef.sigma,
      which           = 2:(K + 1),
      alternative     = "greater",
      R               = 25000,
      mMat            = mMat,
      prob0           = c(Ctrl=1, A=1, B=1, C=1, D=1)[1:(K + 1)],
      stats           = stats.ztest,
      stats.control   = list(sigma = 1),  # known sigma for the z-test
      eff.arm         = eff.arm.simple,
      eff.arm.control = list(lim = mams0$u),
      eff.trial       = eff_trial,
      fut.arm         = fut.arm.simple,
      fut.arm.control = list(lim = mams0$l),
      fut.trial       = fut_trial,
      computation     = "sequential",
      H0              = TRUE,
      extended        = 1,   # store per-trial estimates
      mc.cores        = 1
    )

    # 4. save
    cat("\tSave\n")
    out <- list(mams0 = mams0, mamss0.ztest = mamss0.ztest)
    save(out, file = paste0(path_out, "initial/", id.configw$pos[cw]))

  }
}

cat("\n\t", date(), "\n")
cat("\n\t DONE!\n")
q("no")
