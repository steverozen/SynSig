
set.seed(888)
test.abst <- SignatureAnalyzer4MatchedCatalogs(
  num.runs = 2,
  signatureanalyzer.code.dir =
    "../../data-raw/SignatureAnalzyer.052418/",
  dir.root = "test.random.sigs/",
  slice = 1,
  test.only = TRUE, # Just use a small subset of the input catalog
  overwrite = TRUE,
  mc.cores = 1)
