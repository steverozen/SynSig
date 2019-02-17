## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----readSPprofiles------------------------------------------------------
library(ICAMS)
library(SynSig)

# The column heading for signatures are different in
# SigProfiler output than in PCAWG7, so we turn of
# strict checking.
tmp.read.96 <- function(x) ReadCat96(x, strict = FALSE)

## ------------------------------------------------------------------------
SP.7.analysis <- 
  ReadAndAnalyzeSigs(
    extracted.sigs =
      "res_sp.sp.syn.2018.12.31.v2_signature_patterns_for_7_sigs.csv",
    ground.truth.sigs = "sigProfiler_SBS_signatures_2018_03_28.csv",
    ground.truth.exposures = "sp.syn.exposure.csv",
    read.ground.truth.sigs.fn = tmp.read.96)

SP.7.analysis$avg



## ------------------------------------------------------------------------
knitr::kable(
  SP.7.analysis$match1, 
  caption = 'Best matches from extracted to ground truth'
)

knitr::kable(
  SP.7.analysis$match2,
  caption = 'Best matches from ground truth to extracted'
)


## ---- fig.width = 5------------------------------------------------------

tmp <- 
  lapply(colnames(SP.7.analysis$gt.sigs),
        function(x) {
          PlotCat96(SP.7.analysis$gt.sigs[ ,x, drop = FALSE],
                       type = "signature",
                    grid = FALSE, xlabels = FALSE, cex = 0.6,
                    upper = FALSE)
        })

## ---- fig.width = 5------------------------------------------------------
tmp <- 
  lapply(colnames(SP.7.analysis$ex.sigs),
        function(x) {
          PlotCat96(SP.7.analysis$ex.sigs[ ,x, drop = FALSE],
                       type = "signature",
                    grid = FALSE, xlabels = FALSE, cex = 0.6,
                    upper = FALSE)
        })

## ---- fig.width = 5------------------------------------------------------
SP.9.analysis <- 
  ReadAndAnalyzeSigs(
    extracted.sigs =
      "res_sp.sp.syn.2018.12.31.v2_signature_patterns_for_9_sigs.csv",
    ground.truth.sigs = "sigProfiler_SBS_signatures_2018_03_28.csv",
    ground.truth.exposures = "sp.syn.exposure.csv",
    read.ground.truth.sigs.fn = tmp.read.96)

SP.9.analysis$avg


## ------------------------------------------------------------------------
knitr::kable(
  SP.9.analysis$match1,
  caption = 'Best matches from extracted to ground truth'
)

knitr::kable(
  SP.9.analysis$match2, 
  caption = 'Best matches from ground truth to extracted'
)


## ---- fig.width = 5------------------------------------------------------
tmp <- 
  lapply(colnames(SP.9.analysis$gt.sigs),
        function(x) {
          PlotCat96(SP.9.analysis$gt.sigs[ ,x, drop = FALSE],
                       type = "signature",
                    grid = FALSE, xlabels = FALSE, cex = 0.6,
                    upper = FALSE)
        })

## ---- fig.width = 5------------------------------------------------------
tmp <- 
  lapply(colnames(SP.9.analysis$ex.sigs),
        function(x) {
          PlotCat96(SP.9.analysis$ex.sigs[ ,x, drop = FALSE],
                       type = "signature",
                    grid = FALSE, xlabels = FALSE, cex = 0.6,
                    upper = FALSE)
        })


## ------------------------------------------------------------------------
sessionInfo()

