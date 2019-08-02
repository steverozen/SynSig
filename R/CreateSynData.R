#' @title Extract SynSig parameters for one mutational signature profile
#'
#' @param counts  A vector of mutation counts attributed to one signature across
#'                length(counts) samples. TODO(Steve): rename to exposures
#'
#' @param target.size The length of genomic sequence from which the counts
#'                were derived, in megabases. Deprecated, set this to 1.
#'
#' @return
#'   A 3-element vector with names "prevalence", "mean", and "stdev"
#'
#' @importFrom stats sd
#'
#' @keywords internal

SynSigParamsOneSignature <- function(counts, target.size = 1 ) {

  prevalence <-  length(counts[counts >= 1 ]) / length(counts)

  counts.per.mb <- counts[counts >= 1 ] / target.size

  ## generate log10(mut/mb) values for mean and sd
  mean.per.mb <- mean(log10(counts.per.mb))
  sd.per.mb <- sd(log10(counts.per.mb))

  c(prob=prevalence, mean=mean.per.mb, stdev=sd.per.mb)
}


#' @title Determine 3 parameters for synthetic tumors from an exposure matrix.
#'
#' @param exposures A matrix in which each column is a sample and each row is a mutation
#'         signature, with each element being the "exposure",
#'         i.e. mutation count attributed to a
#'         (sample, signature) pair.
#'
#' @param verbose If > 0 cat various messages.
#'
#' @return A data frame with one row for
#' each of a subset of the input signatures
#'  and the following columns. Signatures not present in
#'  \code{exposures} or present only in a single tumor in
#'  \code{exposures} are removed.
#'
#' \enumerate{
#' \item the proportion of tumors with the signature
#' \item mean(log_10(mutations.per.Mb))
#' \item stdev(log_10(mutations.per.Mb))
#' }
#'
#' @export

GetSynSigParamsFromExposures <-
  function(exposures, # target.size = 1,
           verbose = 0) {
    # target.size Deprecated - was the length of sequence
    # from which the counts were derived; in the future
    # make any adjustments (e.g exome to genome) before
    # providing the exposure matrix.

  integer.counts <- round(exposures, digits=0)
  integer.counts <- integer.counts[rowSums(integer.counts) >0 , ]
  ret1 <- apply(X=integer.counts,
                MARGIN=1,
                FUN=SynSigParamsOneSignature)
                #, target.size)

  # Some standard deviations can be NA (if there is only one tumor
  # with mutations for that signature). We pretend we did not see
  # these signatures. TODO(Steve): impute from similar signatures.
  if(any(is.na(ret1['stdev', ]))) {
    if (verbose > 0) {
      cat("\nWarning, some signatures present in only one sample, dropping:\n")
      cat(colnames(ret1)[is.na(ret1['stdev', ])], "\n")
    }
  }
  retval <- ret1[,!is.na(ret1['stdev',]) , drop = FALSE]
  if (ncol(retval) == 0) {
    stop("No signatures with usable parameters (> 1 sample with exposure)")
  }
  return(retval)
}

#' @title Write SynSig parameters --prevalence, mean(log(exposure))
#'  and sd(log(exposure)) to a file.
#'
#' @param params The parameters to write.
#'
#' @param file   The path to the file to write.
#'
#' @param append Whether to append to or overwrite \code{file} if it already
#'        exists.
#'
#' @param col.names If NA, add column names.
#'
#' @importFrom utils write.table
#'
#' @export

# Needs to be exported for e.g. Create.pancreas.Rmd

WriteSynSigParams <-
  function(params, file, append = FALSE,
           col.names = ifelse(append, FALSE, NA)) {
    # Suppress warning about writing column names
    # on an append.
    suppressWarnings(
      write.table(x = as.data.frame(params),
                  file = file,
                  sep = ",",
                  col.names = col.names,
                  row.names = TRUE,
                  append = append)
    )
  }

#' @title Create synthetic exposures based given parameters
#'
#' @return A matrix with the rows being each signature and the columns being
#' generated samples. Each entry is the count of mutations due to one
#' signature in one sample.
#'
#' @param sig.params Parameters from \code{synsig.parameters.from.attribution}
#'
#' @param num.samples Number of samples to generate
#'
#' @param name Prefix for sample identifiers in the simulated dataset
#'
#' @export

GenerateSyntheticExposures <-
  function(sig.params,
           num.samples = 10,
           name = 'synthetic') {

    sigs <- colnames(sig.params)
    stopifnot(!is.null(sigs))
    prev.present <- sig.params['prob', ] # Note, get a vector
    sig.burden <- sig.params['mean', , drop = F]
    sig.sd <- sig.params['stdev', , drop = F]

    sig.present <- present.sigs(num.samples, prev.present)

    colnames(sig.present) <- paste(name, seq(1, num.samples), sep='.')

    # Create a synthetic exposures for each column (sample)
    # in sig.present.
    retval <-
      apply(sig.present,
            2,
            GenerateSynExposureOneSample,
            sigs,
            sig.burden, ## burden is in mutation per megabase
            sig.sd)
    return(retval)
  }


#' Randomly assign present / absent to each of a set of signatures, and
#' keep trying until at least one is present
#'
#' @param prev.present  Vector of prevalences,
#'   each the prevalence of 1 mutational
#'   signature.
#'
#' @param verbose If > 0, cat some possibly informative messages
#'
#' @return a vector of 0s and 1s of length
#' \code{length(prev.present)}, and for which
#' \code{sum(prev.present) > 0}.
#'
#' @keywords internal
AssignPresentAbsentOneSample <- function(prev.present, verbose = 0) {
  v <- numeric(length(prev.present))
  while (sum(v) < 1) {
    v <- rbinom(length(prev.present), 1, prev.present)
  }
  if (verbose > 0)
    cat("\nAssignPresentAbsentOneSample returning ", v, "\n\n")
  return(v)
}


#' @title Decide which signatures will be present in
#'  the catalogs of synthetic tumors.
#'
#' @param num.tumors Number of tumors to generate
#'
#' @param prev.present Vector of prevalences,
#'   each the prevalence of 1 mutational
#'   signature. The names are the names of the
#'   signatures.
#'
#' @return A matrix with one row
#' per signature and one column per tumor, with
#' 1 in a cell indicated the presence of a signature
#' and 0 indicating absence.
#'
#' @keywords internal

present.sigs <-
  function(num.tumors, prev.present){

    num.process <- length(prev.present)

    present.list <- list()

    for (i in 1:num.process){
      present.each <- rbinom(num.tumors,
                             1,
                             prob = prev.present[i])
      present.list[[i]] <-  present.each
    }

    present <- do.call(rbind,present.list)

    rownames(present) <- names(prev.present)

    # If the column for one tumor has only 0s,
    # re-sample until there is at least on non-0
    # signature.
    for (tumor in 1:ncol(present)){
      if (all(present[,tumor] == rep(0,length(nrow(present))))){
        # present['SBS1',tumor] = 1
        present[ , tumor] <-
          AssignPresentAbsentOneSample(prev.present)
      }
    }
    return(present)
  }


#' @title Using parameters given to generate exposures for synthetic tumors
#'
#' @param tumor Signature presence matrix or exposure matrix for a tumor.
#' It has only one row, and K (# of signatures) columns.
#' Value in each column is the presence flag for a mutational signature:
#' the value can be non-zero(signature is present) or 0(absent).
#' The name of each column should be the name of a signature.
#'
#' @param sig.interest Names of mutational signatures you want to use to
#' generate exposures. It can be all, or part of signatures in colnames(tumor).
#'
#' @param burden.per.sig Mean mutation burden in log10(muts/mb).
#' (counts of mutations per megabase of a tumor sequence)
#' It has one row, and K columns. Each column name refers to a mutational signature.
#'
#' @param sd.per.sig standard deviation of mutation burden.
#' It has one row, and K columns. Each column name refers to a mutational signature.
#'
#' @details ??Determine the intensity of each
#' mutational signature in a tumor, returning mutations per mb
#' using the mean mutation burden per signature and the std dev
#'
#' @importFrom stats rbinom rnorm
#'
#' @keywords internal
GenerateSynExposureOneSample <-
  function(tumor,
           sig.interest,
           burden.per.sig,
           sd.per.sig
  ) {

    ## starts with individual tumors, only generate exposures for signatures
    ## with a flag does not equal to 0.
    active.sigs <- which(tumor != 0)

    for (sigs in active.sigs) {
      stdev <- sd.per.sig[,sigs]
      burden <- burden.per.sig[,sigs]

      ## if std dev is too big, >= 3, max = 3
      ### consider handling this different. the worry is that the variation
      ##  is too large, the sampled mutation burden will be very high,
      ### which will have a mutation burden that is not biologically possible
      if (stdev >= 3) {
        cat("Very large stdev", stdev, "\n")
        stdev = 3
      }

      ## mutational intensity follows a log normal distibution
      ## use the normal distribution with log-ed values instead
      mutations.per.mb <- 10^(rnorm(1,
                                    sd = stdev,
                                    mean = burden))

      tumor[sigs] <- mutations.per.mb
    }

    tumor <- as.matrix(tumor)
    names(tumor) <- sig.interest
    return(tumor)

  }

#' @title Generate synthetic spectra catalogs given
#' signature profiles and synthetic exposures.
#'
#' @param signatures The signature profiles.
#'
#' @param exposures The synthetic exposures.
#'
#' @param sample.id.suffix A string for adding a suffix to
#'  sample ID. For example, if sample.id.suffix is "abc",
#'  then SomeCancerType::s1.33 is changed to
#'  SomeCancerType::s1-abc.33. Actually, this just replaces
#'  the first "." in the sample id with "-" concatenated
#'  to sample.id.suffix. TODO(Steve): probably drop this
#'
#' @return A list of three elements that comprise the
#' synthetic data: \enumerate{
#'  \item \code{ground.truth.catalog}: Spectra catalog for
#'  the software input.
#'  \item \code{ground.truth.signatures}: Signatures active
#'  in \code{ground.truth.catalog}.
#'  \item \code{ground.truth.exposures}: Exposures of \code{ground.truth.signatures}
#'  in \code{ground.truth.catalog}.
#' }
#'
#' @export

CreateSynCatalogs <-
  function(signatures, exposures, sample.id.suffix = NULL) {

  if (any(colSums(exposures) < 1)) warning("Some exposures < 1")

  exposed.sigs <- rownames(exposures)

  # It is an error if there are signatures in exposures that are not
  # in signatures.
  stopifnot(setequal(setdiff(exposed.sigs, colnames(signatures)), c()))

  # VERY IMPORTANT, the next statement guarantees that
  # the order of signatures in rows of exposures is the same as
  # the order of columns in signatures. In addition,
  # it ensure that signatures contains only signatures
  # that are present in exposures.
  #
  signatures <- signatures[ , exposed.sigs]

  # TODO(Steve): in a future versions, for each
  # synthetic tumor, for each signature, accounting
  # for e mutations in that tumor, sample e mutations
  # from the signature profile.
  catalog <- signatures %*% exposures
  stopifnot(!any(colSums(catalog) < 1))

  i.cat <- round(catalog, digits = 0)

  # The condition below can occur even if any(colSUms(catalog) < 1) is FALSE,
  # if before rounding mutiple mutational classes had < 0.5 mutations.
  if (any(colSums(i.cat) == 0))
    warning("Some tumors with 0 mutations")

  if (!is.null(sample.id.suffix)) {
    newcolnames <-
      gsub(".", paste0("-", sample.id.suffix, "."),
           colnames(i.cat), fixed = TRUE)
    stopifnot(colnames(i.cat == colnames(exposures))) #NEW
    colnames(i.cat) <- newcolnames
    colnames(exposures) <- newcolnames
  }

  ## Convert ground.truth.catalog and ground.truth.signatures
  ## into ICAMS acceptable catalogs before outputting the list
  i.cat <- ICAMS::as.catalog(object = i.cat,
                             ref.genome = "hg19",
                             region = "genome",
                             catalog.type = "counts")
  signatures <- ICAMS::as.catalog(object = signatures,
                                  ref.genome = "hg19",
                                  region = "genome",
                                  catalog.type = "counts.signature")

  ## Return a list with:
  ## $ground.truth.catalog: Spectra catalog for the software input
  ## $ground.truth.signatures: Signatures active in $ground.truth.catalog
  ## $ground.truth.exposures: Exposures of $ground.truth.signatures in
  ## $ground.truth.catalog.
  return(list(ground.truth.catalog=i.cat,
              ground.truth.signatures=signatures,
              ground.truth.exposures=exposures))
  # TODO(Steve) In future, add noise
}

#' Merge 2 exposure matrices
#'
#' @param exp1 An exposure matrix
#'
#' @param exp2 An exposure matrix
#'
#' @return The column-wise merge of the two input matrices as
#' with all rownames from either matrix preserved and
#' corresponding entries filled with 0s.
#'
#' @keywords internal
Merge2Exposures <- function(exp1, exp2) {
  # Rows are signatures
  exp.m <- merge(exp1, exp2, by = 0, all = TRUE)
  rownames(exp.m) <- exp.m[ ,1]
  exp.m <- exp.m[ , -1]
  exp.m[is.na(exp.m)] <- 0
  return(as.matrix(exp.m))
}

#' Merge all exposure matrices in a list of matrices
#'
#' @param list.of.exposures A list of exposure matrices
#'
#' @return The column-wise merge of all the input matrices
#' with all rownames from all matrices preserved and
#' corresponding entries filled with 0s.
#'
#' @export

MergeExposures <- function(list.of.exposures) {
  stopifnot(length(list.of.exposures) > 0)
  if (length(list.of.exposures) == 1) {
    return(as.matrix(list.of.exposures[[1]]))
  }
  start <- list.of.exposures[[1]]
  for (i in 2:length(list.of.exposures)) {
    start <- Merge2Exposures(start, list.of.exposures[[i]])
  }
  start2 <- start[order(rownames(start)), ]
  return(as.matrix(start2))
}

#' Create file names in a given directory
#'
#' The directory is provided by the global
#' variable \code{OutDir.dir},
#' which \strong{must} be set by the user. If
#' \code{OutDir.dir} is NULL then just return
#' \code{file.name}.
#'
#' @param file.name The name of the that will be
#' prefixed by \code{OutDir.dir}.
#'
#' @return \code{file.name} prefixed by \code{OutDir.dir}.
#'
#' @export
#'
OutDir <- function(file.name) {
  if (is.null(OutDir.dir)) return(file.name)
  tmp <- OutDir.dir
  n <- nchar(tmp)
  if (substr(tmp, n, n) != "/")
    tmp <- paste0(tmp, "/")
  if (!dir.exists(tmp)) {
    dir.create(tmp)
  }
  return(paste0(tmp, file.name))
}

#' Globally set the location used by \code{OutDir}.
#'
#' @param dir The location of the top level directory
#'
#' @param overwrite If TRUE allow overwriting of existing directory.
#'
#' @param recursive Allow creating directories recursively.
#'
#' @export
SetNewOutDir <- function(dir,
                         overwrite = FALSE,
                         recursive = TRUE) {
  if (dir.exists(dir)) {
    if (overwrite) {
      warning("\nOverwriting ", dir)
    } else stop(dir, " already exists")
  } else {
    dir.create(dir, recursive = recursive)
  }

  ## To assign globally, don't use `<<-`.
  ## Otherwise it won't pass the devtools::check().
  assign("OutDir.dir", "dir", envir = .GlobalEnv)
}

#' Generate synthetic exposures from abstract parameters.
#'
#' Checkpoints the parameters and the synthetic
#' exposures to files. It also checks that the parameters
#' inferred from the synthetic data approximate those
#' inferred from \code{real.exp}.
#'
#' @param parms The actually exposures upon which to base
#' the parameters and synthetic exposures.
#'
#' @param num.syn.tumors Generate this number of synthetic tumors.
#'
#' @param file.prefix Prepend this to output filenames
#'  to indicate the organization of the data.
#'
#' @param sample.id.prefix Prefix for sample identifiers for the
#' synthetic samples.
#'
#' @return A list with elements:
#' \enumerate{
#'  \item \code{parms} The parameters inferred from \code{real.exp}.
#'  \item \code{syn.exp} The synthetic exposures generated from \code{parms}.
#' }
#'
#' @keywords internal

GenerateSynAbstract <-
  function(parms, num.syn.tumors, file.prefix, sample.id.prefix) {
    stopifnot(!is.null(parms))
    froot <- OutDir(file.prefix)

    parm.file <- paste0(froot, ".parms.csv")
    cat("# Original paramaters\n", file = parm.file)
    suppressWarnings( # Suppress warning on column names on append
      WriteSynSigParams(parms, parm.file, append = TRUE,
                      col.names = NA))

    syn.exp <-
      GenerateSyntheticExposures(parms, num.syn.tumors, sample.id.prefix)

    WriteExposure(syn.exp, paste0(froot, ".generic.syn.exposure.csv"))

    # Sanity check
    check.params <- GetSynSigParamsFromExposures(syn.exp)

    # sa.check.param should be similar to parms
    cat("# Parameters derived from synthetic exposures\n",
        file = parm.file, append = TRUE)
    suppressWarnings(
      WriteSynSigParams(check.params, parm.file, append = TRUE))

    missing.sig.names <- setdiff(colnames(parms), colnames(check.params))
    if (length(missing.sig.names) > 0) {
      cat("# Some signatures not represented in the synthetic data:\n",
          file = parm.file, append =  TRUE)
      cat("#", missing.sig.names, "\n", file = parm.file, append = TRUE)
      check.param2 <- matrix(NA, nrow=dim(parms)[1], ncol = dim(parms)[2])
      dimnames(check.param2) <- dimnames(parms)
      check.param2[ , colnames(check.params)] <- check.params
      check.params <- check.param2
    }

    cat("# Difference between original parameters and parameters",
        "derived from synthetic exposures\n",
        file = parm.file, append = TRUE)
    WriteSynSigParams(parms - check.params, parm.file,
                      append = TRUE)
    cat("# The difference should be small\n",
        file = parm.file, append = TRUE)

    return(list(parms=parms, syn.exp=syn.exp))
  }

#' Generate synthetic exposures from real exposures.
#'
#' Checkpoints the parameters and the synthetic
#' exposures to files. It also checks that the parameters
#' inferred from the synthetic data approximate those
#' inferred from \code{real.exp}.
#'
#' @param real.exp The actual (real) exposures upon which to base
#' the parameters and synthetic exposures.
#'
#' @param num.syn.tumors Generate this number of synthetic tumors.
#'
#' @param file.prefix Prepend this to output filenames
#'  to indicate the organization of the data.
#'
#' @param sample.id.prefix Prefix for sample identifiers for the
#' synthetic samples.
#'
#' @return A list with elements:
#' \enumerate{
#'  \item \code{parms} The parameters inferred from \code{real.exp}.
#'  \item \code{syn.exp} The synthetic exposures generated from \code{parms}.
#' }
#'
#' @export
GenerateSynFromReal <-
  function(real.exp, num.syn.tumors, file.prefix, sample.id.prefix) {

  parms <- GetSynSigParamsFromExposures(real.exp)

  WriteExposure(real.exp,
                paste0(OutDir(file.prefix), ".real.input.exposure.csv"))

  return(
    GenerateSynAbstract(
      parms, num.syn.tumors, file.prefix, sample.id.prefix))
  }

#' Create and write a mutational spectra catalog
#'
#' @export
#'
#' @param sigs Signatures to use.
#'
#' @param exp (Synthetic) exposures.
#'
#' @param dir Directory in which to put the signatures;
#' NOTE: this will be a subdirectory based on \code{\link{OutDir}}.
#'
#' @param write.cat.fn Function to write catalogs \strong{or}
#' spectra to files.
#'
#' @param extra.file.suffix Extra string to put before ".csv".
#'
#' @param overwrite If TRUE, overwrite existing directory; useful for
#' debugging / testing.
#'
#' @return Invisibly, the generated catalog.
#'
#' @details Create a file with the catalog \code{syn.data.csv}
#'  and writes \code{sigs} to \code{input.sigs.csv}.
#'
CreateAndWriteCatalog <-
  function(sigs, exp, dir, write.cat.fn, extra.file.suffix = "",
           overwrite = FALSE) {
    info <- CreateSynCatalogs(sigs, exp)

    if (dir.exists(OutDir(dir))) {
      if (!overwrite) stop("\nDirectory ", OutDir(dir), " exists\n")
    } else {
      dir.create(OutDir(dir))
    }

    if (extra.file.suffix == "") {
      suffix <- ".csv"
    } else {
      suffix = paste0(".", extra.file.suffix, ".csv")
    }
    write.cat.fn(info$ground.truth.signatures,
                  OutDir(paste0(dir, "/ground.truth.syn.sigs", suffix)))

    zero.mutation <- which(colSums(info$ground.truth.catalog) == 0)

    if (length(zero.mutation) > 0) {
      warning("Tumors with no mutation:\n\n",
              colnames(info$ground.truth.catalog)[zero.mutation],
              "in", OutDir(dir))
    }
    write.cat.fn(info$ground.truth.catalog,
                 OutDir(paste0(dir, "/ground.truth.syn.catalog", suffix)))
    WriteExposure(info$ground.truth.exposures,
                  OutDir(paste0(dir, "/ground.truth.syn.exposures", suffix)))
    invisible(info$ground.truth.catalog)
  }

#' @keywords internal
MustCreateDir <- function(dir) {
  if (!dir.create(dir, recursive = TRUE)) {
    stop("Unable to create dir ", dir )
  }
}

#' Create and write a mutational spectra catalog
#'
#' @export
#'
#' @param sigs Signatures to use.
#'
#' @param exp (Synthetic) exposures.
#'
#' @param dir Directory in which to put the signatures;
#' NOTE: this will be a subdirectory based on \code{\link{OutDir}}.
#'
#' @param extra.file.suffix Extra string to put before ".csv".
#'
#' @param overwrite If TRUE, overwrite existing directory; useful for
#' debugging / testing.
#'
#' @return Invisibly, the generated catalog.
#'
#' @details Create a file with the catalog \code{syn.data.csv}
#'  and writes \code{sigs} to \code{input.sigs.csv}.
#'
NewCreateAndWriteCatalog <-
  function(sigs, exp, dir, extra.file.suffix = "",
           overwrite = FALSE) {
    info <- CreateSynCatalogs(sigs, exp)

    if (dir.exists(dir)) {
    if (!overwrite) stop("\nDirectory ", dir, " exists\n")
    } else {
      MustCreateDir(dir)
    }

    if (extra.file.suffix == "") {
      suffix <- ".csv"
    } else {
      suffix = paste0(".", extra.file.suffix, ".csv")
    }
    ICAMS::WriteCatalog(info$ground.truth.signatures,
                        paste0(dir, "/ground.truth.syn.sigs", suffix))

    zero.mutation <- which(colSums(info$ground.truth.catalog) == 0)

    if (length(zero.mutation) > 0) {
      warning("Tumors with no mutation:\n\n",
              colnames(info$ground.truth.catalog)[zero.mutation],
              "in", dir)
    }
    ICAMS::WriteCatalog(info$ground.truth.catalog,
                        paste0(dir, "/ground.truth.syn.catalog", suffix))
    WriteExposure(info$ground.truth.exposures,
                  paste0(dir, "/ground.truth.syn.exposures", suffix))
    invisible(info$ground.truth.catalog)
  }

#' Add an R script to a particular subdirectory.
#'
#' @param maxK The \code{maxK} argument for SignatureAnalyzers.
#'
#' @param slice Which subdirectory to put the script into.
#'
#' @param dir.name The name of the subdirectory.
#'
#' @keywords internal

AddScript <- function(maxK, slice,
                      dir.name = c("sa.sa.96",
                                   "sp.sp",
                                   "sa.sa.COMPOSITE",
                                   "sp.sa.COMPOSITE")) {

  lines <- c(
    "",
    "",
    "",
    "library(SynSig)",
    "library(ICAMS)",
    "cat(\"\\n\\nRunning, maxK.for.SA is\", maxK.for.SA, \"\\n\\n\")",
    "RNGkind(kind = \"L'Ecuyer-CMRG\")",
    "set.seed(888)",
    "",
    "reval <- SignatureAnalyzer4MatchedCatalogs(",
    "  num.runs = 20,",
    "  signatureanalyzer.code.dir = \"/home/gmssgr/bin/SignatureAnalzyer.052418/\",",
    "  dir.root = \"..\",",
    "",
    "  overwrite = FALSE,",
    "  maxK = maxK.for.SA,",
    "  mc.cores = 20",
    "  )"
  )

  out.script.name <- paste0(slice, ".run.SA.R")
  lines[1]  <-
    paste0("# Put this file in <top.level.dir>/", dir.name,
           " and run Rscript ", out.script.name)
  lines[2]  <- paste0("maxK.for.SA <- ", maxK)
  lines[14] <- paste0("  slice = ", slice, ",")
  out.name  <- OutDir(paste0(dir.name, "/", out.script.name))
  writeLines(lines, con = out.name)
}


#' Create scripts to run SignatureAnalyzer in all subdirectories
#' of \code{OutDir()}.
#'
#' @param maxK The \code{maxK} argument for SignatureAnalyzer.
#'
#' @export

AddAllScripts <- function(maxK = 30) {
  AddScript(maxK = maxK, 1, "sa.sa.96")
  AddScript(maxK = maxK, 2, "sp.sp")
  AddScript(maxK = maxK, 3, "sa.sa.COMPOSITE")
  AddScript(maxK = maxK, 4, "sp.sa.COMPOSITE")
}


