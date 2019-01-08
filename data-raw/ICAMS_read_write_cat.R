# Functions to read and write catalogs from and to files

ReadCat96 <- function(path, strict = TRUE) {
  # Read 96-channel spectrum or signatures in PCAWG7 format.
  # 
  # Args:
  #   path:   Path to a catalog on disk in the "PCAWG7" format.
  #   strict: If TRUE, do additional checks on the input, and stop
  #           if the checks fail. 
  #   
  # Returns:
  #   A catalog in canonical in-memory format.
  
  cos <- fread(path)
  stopifnot(nrow(cos) == 96)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Mutation type", "Mutation.type"))
    stopifnot(names(cos)[2] == "Trinucleotide")
  }
  ref.gt.var       <- unlist(cos[, 1])
  before.ref.after <- unlist(cos[, 2])
  var <- substring(ref.gt.var, 3, 3)
  out <- cos[ , -(1 : 2)]
  out <- as.matrix(out)
  rownames(out) <- paste0(before.ref.after, var)
  if (strict) {
    stopifnot(rownames(out) == .catalog.row.order96)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[.catalog.row.order96, ]
  return(out)
}

ReadCat192 <- function(path, strict = TRUE) {
  # Read a 192 SNS catalog from path.
  # 
  # Args:
  #   path:   Path to a catalog on disk in the "PCAWG7" format.
  #   strict: If TRUE, do additional checks on the input, and stop
  #           if the checks fail. 
  #   
  # Returns:
  #   A catalog in canonical in-memory format.
  
  cos <- fread(path)
  # cos.copy <- cos # For debugging, testing
  stopifnot(nrow(cos) == 192)
  if (strict) {
    stopifnot(names(cos)[2] %in% c("Mutation type", "Mutation.type"))
    stopifnot(names(cos)[3] == "Trinucleotide")
    stopifnot(names(cos)[1] == "Strand")
  }
  ref.gt.var       <- unlist(cos[, 2])
  before.ref.after <- unlist(cos[, 3])
  
  ## Find the rows labeled with "T", indicating the
  ## SNS is on the transcribed (which is the *antisense*) strand.
  transcribed.strand.pos <- which(cos[, 1] == 'T')
  
  before.ref.after[transcribed.strand.pos] <-
    revc(before.ref.after[transcribed.strand.pos])
  
  var <- substring(ref.gt.var, 3, 3)
  var[transcribed.strand.pos] <- revc(var[transcribed.strand.pos])
  
  tmp <- paste0(before.ref.after, var)
  if (strict) {
    stopifnot(tmp == .catalog.row.order192)
  }
  out <- cos[ , -(1 : 3), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- tmp
  out <- out[.catalog.row.order192, ]
  return(out)
}

ReadCat1536 <- function(path, strict = TRUE) {
  # Read a 1536 SNS catalog from path.
  # 
  # Args:
  #   path:   Path to a catalog on disk in the "PCAWG7" format.
  #   strict: If TRUE, do additional checks on the input, and stop
  #           if the checks fail. 
  #   
  # Returns:
  #   A catalog in canonical in-memory format.
  
  cos <- fread(path)
  stopifnot(nrow(cos) == 1536)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Mutation type", "Mutation.type"))
    stopifnot(names(cos)[2] == "Pentanucleotide")
  }
  names(cos)[1:2] <- c("Mutation type", "Pentanucleotide")
  ref.gt.var       <- cos[["Mutation type"]]
  before.ref.after <- cos[["Pentanucleotide"]]
  var <- substring(ref.gt.var, 3, 3)
  out <- as.matrix(cos[ , -(1 : 2)])
  rownames(out) <- paste0(before.ref.after, var)
  if (strict) {
    stopifnot(rownames(out) == .catalog.row.order1536)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[.catalog.row.order1536, ]
  return(out)
}

ReadCatDNS78 <- function(path, strict = TRUE) {
  # Read a 78 DNS catalog from path.
  # 
  # Args:
  #   path:   Path to a catalog on disk in the "PCAWG7" format.
  #   strict: If TRUE, do additional checks on the input, and stop
  #           if the checks fail. 
  #   
  # Returns:
  #   A catalog in canonical in-memory format.
  
  cos <- fread(path)
  stopifnot(nrow(cos) == 78)
  if (strict) {
    stopifnot(names(cos)[1 : 2] == c("Ref", "Var"))
  }
  names(cos)[1 : 2] <- c("Ref", "Var")
  out <- cos[ , -(1 : 2)]
  out <- as.matrix(out)
  rn <- paste0(cos$Ref, cos$Var)
  diff1 <- sort(setdiff(rn, .catalog.row.order.DNS.78))
  if ( (length(diff1) > 0) 
       && 
       (diff1 == c("CGAA", "CGAC", "CGGA", "TAAC", "TAAG", "TACC"))
    &&
    (sort(setdiff(.catalog.row.order.DNS.78, rn) ==
          c("CGGT", "CGTC", "CGTT", "TACT", "TAGG", "TAGT")))
  ) {
    cat("using temporary hack for old DNS canonicalization\n")
        # CGAA -> CGTT
    rn[rn == "CGAA"] <- "CGTT"
    
        # CGAC -> CGGT
    rn[rn == "CGAC"] <- "CGGT"
    
        # CGGA -> CGTC
    rn[rn == "CGGA"] <- "CGTC"
    
        # TAAC -> TAGT
    rn[rn == "TAAC"] <- "TAGT"
    
        # TAAG -> TACT
    rn[rn == "TAAG"] <- "TACT"
    
        # TACC -> TAGG
    rn[rn == "TACC"] <- "TAGG"
  }
  rownames(out) <- rn
    if (strict) {
      stopifnot(rownames(out) == .catalog.row.order.DNS.78)
    }
  out <- out[.catalog.row.order.DNS.78, ]
  return(out)
}

ReadCatDNS144 <- function(path, strict = TRUE) {
  # Read a 144 DNS catalog from path.
  # 
  # Args:
  #   path:   Path to a catalog on disk in the "PCAWG7" format.
  #   strict: If TRUE, do additional checks on the input, and stop
  #           if the checks fail. 
  #   
  # Returns:
  #   A catalog in canonical in-memory format.
  
  cos <- fread(path)
  stopifnot(nrow(cos) == 144)
  if (strict) {
    stopifnot(names(cos)[1 : 2] == c("Ref", "Var"))
  }
  names(cos)[1 : 2] <- c("Ref", "Var")
  out <- cos[ , -(1 : 2)]
  out <- as.matrix(out)
  rn <- paste0(cos$Ref, cos$Var)
  rownames(out) <- rn
  if (strict) {
    stopifnot(rownames(out) == .catalog.row.order.DNS.144)
  }
  out <- out[.catalog.row.order.DNS.144, ]
  return(out)
}

ReadCatQUAD136 <- function(path, strict = TRUE) {
  # Read a 136 Quad catalog from path.
  # 
  # Args:
  #   path:   Path to a catalog on disk in the "PCAWG7" format.
  #   strict: If TRUE, do additional checks on the input, and stop
  #           if the checks fail. 
  #   
  # Returns:
  #   A catalog in canonical in-memory format.
  
  cos <- fread(path)
  stopifnot(nrow(cos) == 136)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Quad", "quad", "QUAD"))
  }
  names(cos)[1] <- "Quad"
  out <- cos[ , -1]
  out <- as.matrix(out)
  rownames(out) <- cos$Quad
  if (strict) {
    stopifnot(rownames(out) == .catalog.row.order.QUAD.136)
  }
  out <- out[.catalog.row.order.QUAD.136, ]
  return(out)
}

ReadCatID <- function(path, strict = TRUE) {
  # Read a IN (insertion/deletion) catalog from path.
  # 
  # Args:
  #   path:   Path to a catalog on disk in the "PCAWG7" format.
  #   strict: If TRUE, do additional checks on the input, and stop
  #           if the checks fail. 
  #   
  # Returns:
  #   A catalog in canonical in-memory format.
  
  cos <- fread(path)
  stopifnot(nrow(cos) == 83)
  cn <- names(cos)
  ex.cn <- c("Type", "Subtype", "Indel_size", "Repeat_MH_size")
  # Repeat_MH_size is the size of repeat OR microhomology (MH)
  if (strict) { for (i in 1:4) { stopifnot(cn[i] == ex.cn[i]) } }
  names(cos)[1:4] <- ex.cn
  rn <- apply(cos[ , 1:4], MARGIN=1, paste, collapse=":")
  # View(data.frame(mini=rn, good=.catalog.row.order.ID))
  out <- as.matrix(cos[ , -(1 : 4)])
  rownames(out) <- rn
  if (strict) {
    stopifnot(rownames(out) == .catalog.row.order.ID)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[.catalog.row.order.ID, ]
  return(out)
}

WriteCat <- function(ct, path, num.row, row.order, row.header, strict) {
  mut.categories <- rownames(ct)
  stopifnot(num.row == nrow(ct))
  if (strict) {
    stopifnot(mut.categories == row.order)
  }
  ct <- ct[row.order, ]
  DT <- as.data.table(ct)
  fwrite(cbind(row.header, DT), file = path)
}

WriteCat96 <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 96, .catalog.row.order96, .ct.96.row.headers, strict)
}

WriteCat192 <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 192, .catalog.row.order192, .ct.192.row.headers, strict)
} 

WriteCat1536 <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 1536, .catalog.row.order1536, .ct.1536.row.headers, strict)
}

WriteCatDNS78 <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 78, .catalog.row.order.DNS.78, .ct.DNS78.row.headers, strict)
}

WriteCatDNS144 <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 144, .catalog.row.order.DNS.144, .ct.DNS144.row.headers, strict)
}

WriteCatQUAD136 <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 136, .catalog.row.order.QUAD.136, .ct.QUAD136.row.headers, strict)
}

WriteCatID <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 83, .catalog.row.order.ID, .ct.ID.row.headers, strict)
}

TestReadWriteCat <- function() {
  
  read.fn <-
    c(ReadCat96,
      ReadCat192,
      ReadCat1536,
      ReadCatDNS78,
      ReadCatDNS144,
      ReadCatQUAD136,
      ReadCatID
      )
  
  write.fn <-
    c(WriteCat96,
      WriteCat192,
      WriteCat1536,
      WriteCatDNS78,
      WriteCatDNS144,
      WriteCatQUAD136,
      WriteCatID)

  fl <-
    c("data/BTSG_WGS_PCAWG.96.csv",
      "data/BTSG_WGS_PCAWG.192.csv",
      "data/BTSG_WGS_PCAWG.1536.csv",
      "data/BTSG_WGS_PCAWG.dinucs.csv",
      "data/HepG2.dinucs.144.csv",
      "data/HepG2.quad.136.csv",
      "data/BTSG_WGS_PCAWG.indels.csv"
    )
  
  Test1Cat <- function(my.read, my.write, my.file) {
    ct1 <- my.read(my.file)
    my.write(ct1, "tmp.ct.txt")
    ct2 <- my.read("tmp.ct.txt")
    stopifnot(ct1 == ct2)  
  }
  
  discard <- mapply(Test1Cat, read.fn, write.fn, fl)
  unlink("tmp.ct.txt")
  cat("ok\n")
}
