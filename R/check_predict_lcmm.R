# Simulated-data check for predict_lcmm (see R/predict_class.R). Run: Rscript R/check_predict_lcmm.R
# Simulated data only. Does lcmm class prediction on NEW children work through latrend?
suppressPackageStartupMessages({library(latrend); library(lcmm); library(mclust)})

# Wrapper: class of NEW children under a latrend lcModelLcmmGMM fit.
predict_lcmm <- function(m, newdata) {
  hl <- m@model
  hl$call[[1]] <- as.name("hlme")                 # latrend stores the function object
  key <- unique(newdata$id)
  nd <- newdata; nd$id <- match(nd$id, key)       # hlme needs a numeric subject
  pc <- predictClass(hl, newdata = nd, subject = "id")
  out <- setNames(pc$class, key[pc$id])
  out[key]                                        # NA for any child hlme dropped
}
set.seed(1)
n <- 900; cls <- sample(1:3, n, TRUE, prob = c(.6, .3, .1))
mu <- rbind(c(10, -0.5, 0), c(25, 1, 0.3), c(50, 3, -0.5))
sim <- do.call(rbind, lapply(seq_len(n), function(i) {
  t <- (0:3) + runif(4, -0.5, 0.5); b0 <- rnorm(1, 0, 3); b1 <- rnorm(1, 0, 0.5)
  data.frame(id = sprintf("s%04d", i), time = t,
             y = mu[cls[i], 1] + b0 + (mu[cls[i], 2] + b1) * t + mu[cls[i], 3] * t^2 + rnorm(4, 0, 3))
}))
truth <- setNames(cls, sprintf("s%04d", seq_len(n)))
A_ids <- names(truth)[1:450]; B_ids <- names(truth)[451:900]
A <- sim[sim$id %in% A_ids, ]; B <- sim[sim$id %in% B_ids, ]

for (gmm in c(FALSE, TRUE)) {
  cat("\n=====", if (gmm) "GMM" else "GBTM", "\n")
  meth <- lcMethodLcmmGMM(fixed = y ~ time + I(time^2), mixture = ~ time + I(time^2),
                          random = if (gmm) ~ time else ~ -1, id = "id", time = "time",
                          nClusters = 3, init = "gridsearch", gridsearch.rep = 5,
                          gridsearch.maxiter = 10, gridsearch.parallel = FALSE, seed = 7L)
  m <- latrend(meth, data = A, verbose = FALSE)
  aA <- setNames(as.integer(trajectoryAssignments(m)), ids(m))

  # Route 1: latrend's own generic
  r1 <- tryCatch({ pp <- suppressWarnings(predictPostprob(m, newdata = B)); cat("predictPostprob uniform?", all(abs(pp - 1/3) < 1e-9), "\n"); "skipped" },
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  # Route 2: lcmm::predictClass on the underlying hlme object
  r2 <- tryCatch({ pc <- predict_lcmm(m, B); data.frame(id = names(pc), class = pc) },
                 error = function(e) paste("ERROR:", conditionMessage(e)))

  # Sanity: predictClass on A's own data must reproduce A's modal classes
  self <- predict_lcmm(m, A)
  cat("self-prediction matches fitted classes:", mean(self == aA[names(self)]), "
")
  if (is.character(r1)) cat("predictPostprob:", r1, "\n") else
    cat("predictPostprob rows:", length(r1), " ARI vs truth:", round(adjustedRandIndex(r1, truth[unique(B$id)]), 3), "\n")
  if (is.character(r2)) cat("predictClass:", r2, "\n") else {
    cat("predictClass rows:", nrow(r2), " ARI vs truth:",
        round(adjustedRandIndex(r2$class, truth[as.character(r2$id)]), 3), "\n")
    if (!is.character(r1)) cat("routes agree:", mean(r1 == r2$class[match(unique(B$id), r2$id)]), "\n")
  }
  # Missing-wave child: drop 1 obs for some B children (native-missing case)
  Bm <- B[-which(B$id %in% B_ids[1:50])[c(TRUE, FALSE, FALSE, FALSE)], ]
  pcm <- predict_lcmm(m, Bm)
  cat("with 3-wave children: non-NA", sum(!is.na(pcm)), "of", length(unique(Bm$id)), "\n")
}
