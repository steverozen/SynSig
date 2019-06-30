# Compare SigProfiler and SignatureAnalyzer extracted signatures

res1 <- MatchSigs2Directions(sp.sigs, sa.96.sigs)

dummy.exposure1 <- matrix(10, ncol=1, nrow=ncol(sp.sigs))
rownames(dummy.exposure1) <- colnames(sp.sigs)

res2 <- MatchSigsAndRelabel(
  ex.sigs = sa.96.sigs, gt.sigs = sp.sigs, exposure = dummy.exposure1)

dummy.exposure2 <- matrix(10, ncol=1, nrow=ncol(sa.96.sigs))
rownames(dummy.exposure2) <- colnames(sa.96.sigs)

res3 <- MatchSigsAndRelabel(
  gt.sigs = sa.96.sigs, ex.sigs = sp.sigs, exposure = dummy.exposure2)

file.exists(devtools::package_file("data-raw/compare.SP.and.SA.sigs"))

sa.96.total.exp <- matrix(rowSums(sa.all.real.exposures), ncol = 1)
rownames(sa.96.total.exp) <- rownames(sa.all.real.exposures)
sa.mat <- cbind(res2$match1, sa.96.total.exp[rownames(res2$match1), ])
write.csv(
  sa.mat,
  devtools::package_file(
    "data-raw/compare.SP.and.SA.sigs/SA.to.SP.csv"))

sp.total.exp <- matrix(rowSums(sp.all.real.exposures), ncol = 1)
rownames(sp.total.exp) <- rownames(sp.all.real.exposures)
sp.mat <- cbind(res2$match2, sp.total.exp[rownames(res2$match2), ])
write.csv(
  sp.mat, # res2$match2,
  devtools::package_file(
    "data-raw/compare.SP.and.SA.sigs/SP.to.SA.csv"))

PlotCatSNS96ToPdf(
  catalog = res2$ex.sigs, 
  type = "signature",
  filename =   devtools::package_file(
    "data-raw/compare.SP.and.SA.sigs/SA.to.SP.pdf"))

PlotCatSNS96ToPdf(
  catalog = res3$ex.sigs, 
  type = "signature",
  filename =   devtools::package_file(
    "data-raw/compare.SP.and.SA.sigs/SP.to.SA.pdf"))

misc.file <- devtools::package_file(
  "data-raw/compare.SP.and.SA.sigs/other.info.txt")

cat("average =", res2$avg, "\n", file = misc.file)
cat("no best match1\n", res2$ground.truth.with.no.best.match, "\n", 
    file = misc.file, append = TRUE)
cat("no best match2\n", res3$ground.truth.with.no.best.match, "\n",
    file = misc.file, append = TRUE)






