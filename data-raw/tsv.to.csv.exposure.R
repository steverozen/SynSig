.fix1 <- function(x) {
infile <- paste0(x, "/sp.sp/ground.truth.syn.exposures.txt")
outfile <- paste0(x, "/sp.sp/ground.truth.syn.exposures.csv")
t1 <- read.table(infile, sep = "\t")
write.csv(t1, outfile )
}
lapply(dir(pattern = "Rsq"), .fix1)

library(R.utils)
.fix2 <- function(x) {
  infile <- paste0(x, "/sp.sp/ground.truth.syn.signatures.csv")
  outfile <- paste0(x, "/sp.sp/ground.truth.syn.sigs.csv")
  if (!file.exists(outfile)) {
    copyFile(srcPathname = infile,  destPathname = outfile)
  }
}
lapply(dir(pattern = "Rsq"), .fix2)
