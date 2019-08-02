
## Load required packages
library(ICAMS)
library(SynSig)



## Specify slopes and Rsqs for the datasets
slopes <- c(0.1,0.5,1)
Rsqs <- c(0.1,0.2,0.3,0.6)
datasetNames <- c()
for(slope in slopes)
  for(Rsq in Rsqs)
    datasetNames <- c(datasetNames, paste0("S.",slope,".Rsq.",Rsq))


## Run Extraction and attribution packages
## sigproextractor (Python package) and MultiModalMuSig (Julia package)
## needs to be run with external script.
extrAttrPackages <- c("signeR","sigfit","MutationalPatterns",
                      "SomaticSignatures","hdp")
for(extrAttrPackage in extrAttrPackages){
  func <- get(paste0("Run",extrAttrPackage))
  for(datasetName in datasetNames){
    func(input.catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
         read.catalog.function = ICAMS::ReadCatalog,
         write.catalog.function = ICAMS::WriteCatalog,
         out.dir = paste0(datasetName,"/sp.sp/ExtrAttr/",extrAttrPackage,".results"),
         K.range = c(1,10),
         overwrite = T)
  }
}



## Run Attribution-only packages (deconstructSigs, YAPSA, SignatureEstimation)
## Note: SignatureEstimation has two methods:
attrOnlyPackages <- c("deconstructSigs","YAPSA")

for(attrOnlyPackage in attrOnlyPackages){
  func <- get(paste0("Run",attrOnlyPackage,"AttributeOnly"))
  for(datasetName in datasetNames){
    func(input.catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
         gt.sigs.file = paste0(datasetName,"/sp.sp/ground.truth.syn.sigs.csv"),
         read.catalog.function = ICAMS::ReadCatalog,
         out.dir = paste0(datasetName,"/sp.sp/Attr/",attrOnlyPackage,".results"),
         overwrite = T)
  }
}
