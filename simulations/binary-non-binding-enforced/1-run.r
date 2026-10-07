# Simulate every config of the grid: MAMS design, MAMSS evaluation.
# Usage: Rscript 1-run.r 1 in out
#
# Non-binding boundaries (binding = FALSE); futility is applied in the
# simulation (fut.arm.simple), giving realistic power and ESS.
# Binary endpoint: Bernoulli outcomes (rBI), pooled two-proportion test (stats.ptest).

args <- commandArgs(trailingOnly = TRUE)

library(MAMS)
library(MAMSS)
library(mvtnorm)
library(dotfunctions)
source("../binary-design-helper.R")
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

    # true response rates
    pc <- id.configw$pc[cw]
    pt <- pc + (2 * id.configw$p[cw] - 1)
    if (pt <= 0 || pt >= 1) {
      stop(sprintf("config %s: pi_t = %.3f out of (0,1); adjust pc in 0-par.r",
                   id.configw$pos[cw], pt))
    }
    delta_bin <- pt - pc                       # binary effect (risk difference)

    cat("\tDesign (binary_design(): MAMS boundaries + H1-correct n search, non-binding)\n")

    # 1. design: boundaries from MAMS, sample size from binary_design()
    #    (accounts for the H1 variance of the pooled test)
    mams0 <- binary_design(
      K      = K,
      J      = J,
      alpha  = id.configw$alpha[cw],
      power  = id.configw$power[cw],
      delta  = delta_bin,
      pc     = pc,
      pt     = pt,
      method = id.configw$method[cw],
      r      = 1:J,
      r0     = 1:J,
      ushape     = "obf",
      lshape     = "triangular",
      binding    = FALSE
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

    # 3. operating characteristics with MAMSS (logit scale: control and
    #    null arms at pi_c, T1 at pi_t)
    coef.mu <- c(qlogis(pc), qlogis(pt) - qlogis(pc), rep(0, K - 1))

    # sep: arms stop individually, the trial never stops as a whole
    if (id.configw$method[cw] == "sep") {
      eff_trial <- function(eff.target) FALSE
      fut_trial <- function(fut.target) FALSE
    } else {
      eff_trial <- eff.trial.any
      fut_trial <- fut.trial.all
    }

    cat("\tEvaluation (MAMSS::mamss - rBI, stats.ptest, fut.arm.simple [enforced])\n")

    mamss0.ptest <- mamss(
      model           = outcome ~ intervention,
      var             = list(outcome   = rBI,
                             intervention = alloc.balanced),
      var.control     = list(outcome = list(bd = 1)),  # bd = 1 => Bernoulli
      coef.mu         = coef.mu,
      which           = 2:(K + 1),
      alternative     = "greater",
      R               = 25000,
      mMat            = mMat,
      prob0           = c(Ctrl=1, A=1, B=1, C=1, D=1)[1:(K + 1)],
      stats           = stats.ptest,
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
    out <- list(mams0 = mams0, mamss0.ptest = mamss0.ptest,
                rates = c(pc = pc, pt = pt))
    save(out, file = paste0(path_out, "initial/", id.configw$pos[cw]))

  }
}

cat("\n\t", date(), "\n")
cat("\n\t DONE!\n")
q("no")
