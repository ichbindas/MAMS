# Sample size for binary-endpoint MAMS designs tested with the pooled
# two-proportion test (MAMSS::stats.ptest).
#
# Under H1 the pooled statistic Z has (1:1 allocation, pbar = (pc + pt) / 2)
#   mean(Z) = delta * sqrt(n / (2 * pbar * (1 - pbar)))
#   var(Z)  = R = (pt * (1 - pt) + pc * (1 - pc)) / (2 * pbar * (1 - pbar)) <= 1
# MAMS::mams() assumes var(Z) = 1, which a single 'sd' cannot correct. So the
# boundaries (u, l) and rMat come from mams(sample.size = FALSE), and n is
# searched here with the correct mean and variance.
#
# The search uses the "sep" boundaries, where power reduces to a single-arm
# group-sequential calculation. The same n is used for "simultaneous"; for
# the grid used here the two methods' sample sizes differ by less than 0.3%.
# Boundaries passed to the simulation are always those of the requested method.

binary_boundaries <- function(K, J, alpha, power, delta, method,
                               r = 1:J, r0 = 1:J,
                               ushape = "obf", lshape = "triangular",
                               binding = TRUE, seed = 20251123) {
  set.seed(seed)
  d0 <- MAMS::mams(K = K, J = J, alpha = alpha, power = power,
                    delta = delta, delta0 = 0, sd = 1, r = r, r0 = r0,
                    ushape = ushape, lshape = lshape, method = method,
                    binding = binding, sample.size = FALSE,
                    parallel = FALSE, print = FALSE)
  list(u = d0$u, l = d0$l, rMat = d0$rMat)
}

# Sigma_std: J x J correlation matrix of the standardized GSD test statistics
# across looks, for equal-increment information r (used by MAMS' "sep" path;
# depends only on r, not on K, delta or sd).
sigma_std_looks <- function(r, J) {
  bottom <- matrix(r, J, J)
  top <- matrix(rep(r, rep(J, J)), J, J)
  top[upper.tri(top)] <- t(top)[upper.tri(top)]
  bottom[upper.tri(bottom)] <- t(bottom)[upper.tri(bottom)]
  sqrt(top / bottom)
}

# R-corrected sample size search (per-arm "n" in MAMS' rMat*n convention),
# using the "sep" (single-arm GSD) boundaries u/l and the TRUE mean/variance
# of the pooled test statistic under H1 (see header).
binary_design_n <- function(J, alpha, power, delta, pc, pt, r = 1:J,
                             u_sep, l_sep, nmax = 200000) {
  pbar <- (pc + pt) / 2
  sig_pool <- sqrt(pbar * (1 - pbar))
  Rvar <- (pt * (1 - pt) + pc * (1 - pc)) / (2 * pbar * (1 - pbar))
  Sigma_std <- sigma_std_looks(r, J)

  achieved_power <- function(n) {
    ncp <- delta * sqrt((r * n) / (2 * sig_pool^2))
    if (J == 1) {
      return(pnorm(u_sep[1], mean = ncp[1], sd = sqrt(Rvar), lower.tail = FALSE))
    }
    pi <- pnorm(u_sep[1], mean = ncp[1], sd = sqrt(Rvar), lower.tail = FALSE)
    for (j in 2:J) {
      pi <- pi + mvtnorm::pmvnorm(
        lower = c(l_sep[1:(j - 1)], u_sep[j]), upper = c(u_sep[1:(j - 1)], Inf),
        mean = ncp[1:j], sigma = Rvar * Sigma_std[1:j, 1:j]
      )[1]
    }
    pi
  }

  n <- 1
  while (achieved_power(n) < power && n <= nmax) n <- n + 1
  if (n > nmax) stop("binary_design_n: n search exceeded nmax")
  list(n = n, Rvar = Rvar, sig_pool = sig_pool, achieved = achieved_power(n))
}

# Full design for one (K, J, pc, p, method) config: boundaries for the
# REQUESTED method (used for the actual MAMSS simulation), sample size from
# the R-corrected "sep" search (see scope note above).
binary_design <- function(K, J, alpha, power, delta, pc, pt, method,
                           r = 1:J, r0 = 1:J,
                           ushape = "obf", lshape = "triangular",
                           binding = TRUE, seed = 20251123) {
  b_method <- binary_boundaries(K, J, alpha, power, delta, method,
                                 r, r0, ushape, lshape, binding, seed)
  b_sep <- if (identical(method, "sep")) {
    b_method
  } else {
    binary_boundaries(K, J, alpha, power, delta, "sep",
                       r, r0, ushape, lshape, binding, seed)
  }
  sz <- binary_design_n(J, alpha, power, delta, pc, pt, r,
                         u_sep = b_sep$u, l_sep = b_sep$l)
  list(n = sz$n, u = b_method$u, l = b_method$l, rMat = b_method$rMat,
       Rvar = sz$Rvar, sig_pool = sz$sig_pool, achieved = sz$achieved)
}
