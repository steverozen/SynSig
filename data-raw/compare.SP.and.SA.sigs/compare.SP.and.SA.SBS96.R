# Compare SigProfiler and SignatureAnalyzer extracted signatures

# res1 <- MatchSigs2Directions(sp.sigs, sa.96.sigs)

dummy.exposure1 <- matrix(10, ncol=1, nrow=ncol(sp.sigs))
rownames(dummy.exposure1) <- colnames(sp.sigs)

res1 <- MatchSigsAndRelabel(
  ex.sigs = sa.96.sigs, gt.sigs = sp.sigs, exposure = dummy.exposure1)

dummy.exposure2 <- matrix(10, ncol=1, nrow=ncol(sa.96.sigs))
rownames(dummy.exposure2) <- colnames(sa.96.sigs)

res2 <- MatchSigsAndRelabel(
  gt.sigs = sa.96.sigs, ex.sigs = sp.sigs, exposure = dummy.exposure2)

file.exists(devtools::package_file("data-raw/compare.SP.and.SA.sigs"))

sa.96.total.exp <- matrix(rowSums(sa.all.real.exposures), ncol = 1)
rownames(sa.96.total.exp) <- rownames(sa.all.real.exposures)
sa.mat <- cbind(res1$match1, sa.96.total.exp[rownames(res1$match1), ])
write.csv(
  sa.mat,
  devtools::package_file(
    "data-raw/compare.SP.and.SA.sigs/SA.to.SP.csv"))

sp.total.exp <- matrix(rowSums(sp.all.real.exposures), ncol = 1)
rownames(sp.total.exp) <- rownames(sp.all.real.exposures)
sp.mat <- cbind(res1$match2, sp.total.exp[rownames(res1$match2), ])
write.csv(
  sp.mat, # res1$match2,
  devtools::package_file(
    "data-raw/compare.SP.and.SA.sigs/SP.to.SA.csv"))

PlotCatSNS96ToPdf(
  catalog = res1$ex.sigs,
  type = "signature",
  filename =   devtools::package_file(
    "data-raw/compare.SP.and.SA.sigs/SA.to.SP.pdf"))

PlotCatSNS96ToPdf(
  catalog = res2$ex.sigs,
  type = "signature",
  filename =   devtools::package_file(
    "data-raw/compare.SP.and.SA.sigs/SP.to.SA.pdf"))

misc.file <- devtools::package_file(
  "data-raw/compare.SP.and.SA.sigs/other.info.txt")

cat("average =", res1$avg, "\n", file = misc.file)
cat("no best match1\n", res1$ground.truth.with.no.best.match, "\n",
    file = misc.file, append = TRUE)
cat("no best match2\n", res2$ground.truth.with.no.best.match, "\n",
    file = misc.file, append = TRUE)


# Are SA signatures more sparse?


ineq::Gini(sp.sigs[ , "SBS5"])
# [1] 0.402699
reldist::gini(sp.sigs[ , "SBS5"])
# Same
lawstat::gini.index(sp.sigs[ , "SBS5"])
# Gini Index = 0.40694


ineq::Gini(sa.96.sigs[, "BI_COMPOSITE_SBS5_P"])
# [1] 0.5041189
reldist::gini(sa.96.sigs[, "BI_COMPOSITE_SBS5_P"])
# Same
lawstat::gini.index(sa.96.sigs[, "BI_COMPOSITE_SBS5_P"])
# Gini Index = 0.50943
ineq::Gini(sa.96.sigs[, "BI_COMPOSITE_SBS5_P"], corr = TRUE)
# [1] 0.5094255

ineq::Gini(sp.sigs[ , "SBS3"])
# [1] 0.3244002
ineq::Gini(sa.96.sigs[, "BI_COMPOSITE_SBS3_P"])
# [1] 0.4307233

ineq::Gini(sp.sigs[ , "SBS40"])
# [1] 0.3773224
ineq::Gini(sa.96.sigs[, "BI_COMPOSITE_SBS40_P"])
# [1] 0.4816323

ineq::Gini(sp.sigs[ , "SBS1"])
# [1] 0.9375041
ineq::Gini(sa.96.sigs[, "BI_COMPOSITE_SBS1_P"])
# [1] 0.8242305

ineq::Gini(sp.sigs[ , "SBS7a"])
# [1] 0.9326029

ineq::Gini(sa.96.sigs[, "BI_COMPOSITE_SBS7a_S"])
# [1] 0.9299643

# Another topic -- comparision to COSMICv2

COSMICv2.file <- devtools::package_file(
  "data-raw/COSMICv2_signatures_probabilities.txt")

COSMICv2 <- as.data.frame(data.table::fread(COSMICv2.file))

rownames(COSMICv2) <-
  paste0(COSMICv2[ , 2], substr(COSMICv2[ ,1], 3, 3))

COSMICv2 <- COSMICv2[ , -(1:3)]

COSMICv2 <- COSMICv2[ICAMS::catalog.row.order$SNS96, ]

COSMICv2 <- COSMICv2[ , -(31:37)]

ICAMS::PlotCatSNS96ToPdf(
  COSMICv2,
  devtools::package_file(
    "data-raw/compare.SP.and.SA.sigs/COSMICv2.pdf"),
  type = "signature")

resCOSMIC <- MatchSigsAndRelabel(
  ex.sigs = COSMICv2, gt.sigs = sp.sigs, exposure = dummy.exposure1)

write.csv(
  resCOSMIC$match2,
  devtools::package_file(
    "data-raw/compare.SP.and.SA.sigs/COSMICv3.vs.v2.csv"))
