library(testthat)

nb_design <- function(method, binding, ...) {
  mams(K = 2, J = 2, alpha = 0.05, power = 0.8, r = 1:2, r0 = 1:2,
       p = 0.65, p0 = 0.55, ushape = "obf", lshape = "triangular",
       print = FALSE, method = method, binding = binding, ...)
}

testthat::test_that("binding must be a single logical value", {
  for (m in c("simultaneous", "sep")) {
    testthat::expect_error(
      nb_design(m, binding = "yes", sample.size = FALSE, nsim = NULL),
      "'binding' must be a single logical value"
    )
    testthat::expect_error(
      nb_design(m, binding = c(TRUE, FALSE), sample.size = FALSE, nsim = NULL),
      "'binding' must be a single logical value"
    )
  }
})

testthat::test_that("binding defaults to TRUE and is stored in the result", {
  for (m in c("simultaneous", "sep")) {
    res <- mams(K = 2, J = 2, p = 0.65, p0 = 0.55, sample.size = FALSE,
                print = FALSE, nsim = NULL, method = m)
    testthat::expect_true(res$binding)
    testthat::expect_false(
      nb_design(m, binding = FALSE, sample.size = FALSE, nsim = NULL)$binding
    )
  }
})

testthat::test_that("objects without binding (e.g. from dtl) can be reused", {
  dtl <- mams(K = c(4, 1), J = 2, method = "dtl", print = FALSE, nsim = NULL)
  for (m in c("simultaneous", "sep")) {
    res <- suppressWarnings(mams(dtl, method = m, print = FALSE, nsim = NULL))
    testthat::expect_true(res$binding)
  }
})

testthat::test_that("print marks the lower bound as non-binding", {
  for (m in c("simultaneous", "sep")) {
    out_nb <- capture.output(print(
      nb_design(m, binding = FALSE, sample.size = FALSE, nsim = NULL)))
    out_b <- capture.output(print(
      nb_design(m, binding = TRUE, sample.size = FALSE, nsim = NULL)))
    testthat::expect_true(any(grepl("(non-binding)", out_nb, fixed = TRUE)))
    testthat::expect_false(any(grepl("(non-binding)", out_b, fixed = TRUE)))
  }
})

testthat::test_that("non-binding efficacy boundary is wider than binding", {
  for (m in c("simultaneous", "sep")) {
    set.seed(1)
    b <- nb_design(m, binding = TRUE, sample.size = FALSE, nsim = NULL)
    set.seed(1)
    nb <- nb_design(m, binding = FALSE, sample.size = FALSE, nsim = NULL)
    testthat::expect_gt(nb$u[2], b$u[2])
  }
})

testthat::test_that("non-binding equals binding for a single-stage design", {
  for (m in c("simultaneous", "sep")) {
    set.seed(1)
    b <- mams(K = 2, J = 1, r = 1, r0 = 1, p = 0.65, p0 = 0.55,
              sample.size = FALSE, print = FALSE, nsim = NULL,
              method = m, binding = TRUE)
    set.seed(1)
    nb <- mams(K = 2, J = 1, r = 1, r0 = 1, p = 0.65, p0 = 0.55,
               sample.size = FALSE, print = FALSE, nsim = NULL,
               method = m, binding = FALSE)
    testthat::expect_equal(nb$u, b$u)
  }
})

testthat::test_that("sep simulation removes arms that cross efficacy (non-binding)", {
  set.seed(1)
  res <- nb_design("sep", binding = FALSE, nsim = 1000)
  ess_t1 <- res$sim$H1$main$ess[2, ]
  eff <- res$sim$H1$main$efficacy
  testthat::expect_lt(ess_t1$ess, ess_t1$high)
  testthat::expect_true(all(eff[, ncol(eff)] <= 1))
})
