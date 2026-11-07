# =============================================================================
# R/predict_class.R: assign NEW children to an existing fit's classes.
# Used by the site split-half (prediction strength) analysis.
#   k-means           nearest centre
#   kml, Ward (+ average / complete linkage)
#                     nearest class mean trajectory. Hierarchical clustering has
#                     no native predict; nearest centroid is the standard
#                     approximation and is reported as one.
#   PAM               nearest medoid, in the fit's metric
#   Gaussian mixture  predict.Mclust (MAP class)
#   GBTM, GMM         lcmm posterior class from the training model (predict_lcmm)
# `res` is a fitter's return value (R/fit_task.R); d_train / d_new come from make_d().
# Returned labels are the fit's own labels (before relabelling).
# =============================================================================

nearest <- function(Xn, C, metric = "euclidean") {
  D <- sapply(seq_len(nrow(C)), function(j) {
    dif <- sweep(Xn, 2, C[j, ])
    if (metric == "manhattan") rowSums(abs(dif)) else rowSums(dif^2)
  })
  as.integer(rownames(C))[max.col(-matrix(D, nrow(Xn)), ties.method = "first")]
}

centroids <- function(X, cls) {
  C <- rowsum(X, cls) / as.vector(table(cls))
  rownames(C) <- sort(unique(cls)); C
}

predict_new <- function(method, res, d_train, d_new, args = list()) {
  switch(method,
    kmeans = { C <- res$model; rownames(C) <- seq_len(nrow(C)); nearest(d_new$X, C) },
    kml = , ward = nearest(d_new$X, centroids(d_train$X, res$cls)),
    pam = { C <- res$model; rownames(C) <- seq_len(nrow(C))
            nearest(d_new$X, C, if (is.null(args$metric)) "euclidean" else args$metric) },
    llpa = as.integer(predict(res$model, newdata = d_new$X)$classification),
    gbtm = , gmm = unname(predict_lcmm(res$model, d_new$long)[d_new$ids]),
    stop("no predict for ", method))
}

# GBTM / GMM (latrend lcModelLcmmGMM). Checked on simulated data 2026-10-24
# (self-prediction reproduces the fitted classes exactly; ARI vs truth .98-.99 on
# held-out children, including children with only 3 waves). Three gotchas:
# - latrend's predictPostprob() is NOT implemented for lcModelLcmmGMM. It warns
#   and returns a uniform matrix, so never use it for new data.
# - latrend builds the hlme call with do.call, so call[[1]] is the function object;
#   lcmm::predictClass needs a name there.
# - hlme needs a numeric subject, so ids are mapped to integers and back.
# Also takes an hlme object directly (fit_gmm_conv in R/fit_task.R).
predict_lcmm <- function(m, newdata) {
  hl <- if (inherits(m, "hlme")) m else m@model
  hl$call[[1]] <- as.name("hlme")
  key <- unique(newdata$id)
  nd <- newdata; nd$id <- match(nd$id, key)
  pc <- lcmm::predictClass(hl, newdata = nd, subject = "id")
  out <- setNames(as.integer(pc$class), key[pc$id])
  out[key]                                        # NA for any child hlme dropped
}
