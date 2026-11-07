# =============================================================================
# 00_install.R: one-time package install, with exact versions.
# kml, its dependency longitudinalData, and clv were archived from CRAN on
# 2026-01-16. They are installed from the CRAN archive at the versions pinned below.
# clv 0.3-2.5 includes <R_ext/PrtUtil.h>, which R 4.6 removed. Nothing in clv
# uses it, so the include line is deleted before building (no other changes).
# =============================================================================

cran <- c("latrend", "lcmm", "mclust", "mclustcomp", "cluster", "clusterCrit",
          "mice", "nnet", "ggalluvial", "patchwork", "tidyverse",
          "rgl")                                # rgl: longitudinalData dep
# The CRAN binary of rgl links XQuartz's libGLU, which this Mac lacks, so rgl
# is built from source without OpenGL. kml uses rgl only for 3-D plots.
REPO <- "https://cloud.r-project.org"
if (!"rgl" %in% rownames(installed.packages())) {
  # rgl's own dependencies as binaries first, so only rgl compiles from source.
  db <- available.packages(repos = REPO)
  rgl_deps <- tools::package_dependencies("rgl", db = db, which = c("Depends", "Imports", "LinkingTo"),
                                          recursive = TRUE)[[1]]
  rgl_deps <- intersect(setdiff(rgl_deps, rownames(installed.packages())), rownames(db))
  if (length(rgl_deps)) install.packages(rgl_deps, repos = REPO)
  # --disable-libpng: the second Mac (M5, 2026-10-24) has no libpng headers.
  # kml never uses rgl's PNG snapshots.
  install.packages("rgl", repos = REPO, type = "source", dependencies = FALSE,
                   configure.args = c("--disable-opengl", "--disable-libpng"))
}
need <- setdiff(cran, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = REPO)

archived <- c(clv              = "clv_0.3-2.5.tar.gz",
              longitudinalData = "longitudinalData_2.4.7.tar.gz",
              kml              = "kml_2.5.0.tar.gz")
for (p in names(archived)) {
  if (p %in% rownames(installed.packages())) next
  url <- sprintf("https://cran.r-project.org/src/contrib/Archive/%s/%s", p, archived[[p]])
  tmp <- tempfile(); dir.create(tmp)
  tgz <- file.path(tmp, archived[[p]])
  download.file(url, tgz, quiet = TRUE)
  if (p == "clv") {
    untar(tgz, exdir = tmp)
    h <- file.path(tmp, "clv/src/VCLDefines.h")
    writeLines(grep("R_ext/PrtUtil.h", readLines(h), invert = TRUE, value = TRUE), h)
    tgz <- file.path(tmp, "clv")
  }
  # longitudinalData attaches misc3d only for its 3-D density plots. misc3d
  # loads tcltk, which needs XQuartz, so the dependency is removed from
  # DESCRIPTION and NAMESPACE. kml never calls the 3-D plots.
  if (p == "longitudinalData") {
    untar(tgz, exdir = tmp)
    d <- file.path(tmp, "longitudinalData/DESCRIPTION")
    writeLines(sub(",misc3d", "", readLines(d), fixed = TRUE), d)
    n <- file.path(tmp, "longitudinalData/NAMESPACE")
    writeLines(sub(",misc3d)", ")", readLines(n), fixed = TRUE), n)
    tgz <- file.path(tmp, "longitudinalData")
  }
  install.packages(tgz, repos = NULL, type = "source")
}

ip <- installed.packages()
pk <- c(cran, names(archived))
stopifnot(all(pk %in% rownames(ip)))
writeLines(sprintf("%s %s", pk, ip[pk, "Version"]), "output/logs/package_versions.txt")
