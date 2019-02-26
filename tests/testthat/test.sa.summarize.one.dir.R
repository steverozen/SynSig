

SummarizeSigOneSACOMPOSITESubdir(
  third.level.dir = "./test.sa.sa.COMPOSITE/sa.results/",
  ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
  which.run = "1.run.sa.sa.COMPOSITE")

SummarizeSigOneSA96Subdir(
  third.level.dir = "./test.sa.sa.96/sa.results/",
  ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
  which.run = "1.run.sa.sa.96",
  write.png = FALSE) # ICAMS version skew / Set to false for some OS?
