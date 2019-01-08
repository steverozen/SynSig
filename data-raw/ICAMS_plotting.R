###############################################################################
# Plotting functions for SNS 96, 192 and 1536 catalog start here
###############################################################################

PlotCat96 <- function(catalog, id, type = "density", abundance = NULL) {
  # Plot the SNS 96 mutation catalog of one sample.
  # 
  # Args:
  #   catalog:   A matrix whose rownames indicate the 96 SNS mutation types 
  #              while its columns contain the counts of each mutation type.
  #   id:        The ID information of the sample which has mutations.
  #   type:      A value indicating the type of the graph. 
  #              If type = "density", the graph will plot the rates of 
  #              mutations per million trinucleotides for each mutation type.
  #              If type = "counts", the graph will plot the occurrences of   
  #              the 96 mutation types in the sample.
  #              If type = "signature", the graph will plot mutation 
  #              signatures of the sample.
  #              The default value for type is "density".
  #   abundance: A matrix containing trinucleotide abundance information.
  #              To be used only when type = "density".
  
  stopifnot(dim(catalog) == c(96, 1))
  stopifnot(rownames(catalog) == .catalog.row.order96)
  
  class.col <- c("#0000ff",  # dark blue
                 "#000000",  # black
                 "#ff4040",  # red
                 "#838383",  # grey
                 "#40ff40",  # green
                 "#ff667f")  # pink
  cols <- rep(class.col, each = 16)
  maj.class.names <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  num.classes <- length(catalog)
  
  if (type == "density") {
    # Calculate rate of mutations per million trinucleotides for the catalog
    rate <- double(96)
    for (i in 1 : 96) {
      rate[i] <- 
        catalog[i] * 1000000 / abundance[substr(rownames(catalog)[i], 1, 3), ]
    }
    
    # Get ylim
    ymax <- max(rate)
    
    # Barplot
    bp <- barplot(rate, xaxt = "n", yaxt = "n", xaxs = "i", 
                  xlim = c(-1, 230), lwd = 3, space = 1.35, border = NA, 
                  col = cols, ylab = "mut/million")
    
    # Write the mutation counts on top of graph
    for (i in 1 : 6) {
      j <- 13 + 16 * (i - 1)
      k <- 1 + 16 * (i - 1)
      text(bp[j], ymax * 1.15, labels = sum(catalog[k : (16 * i), ]), 
           xpd = NA, cex = 0.8)
    } 
  } 
  
  if (type == "counts") {
    # Get ylim
    ymax <- max(catalog[, 1])
    
    # Barplot
    bp <- barplot(catalog[, 1], xaxt = "n", yaxt = "n", xlim = c(-1, 230),
                  xaxs = "i", lwd = 3, space = 1.35, border = NA, 
                  col = cols, ylab = "counts")
    
    # Write the mutation counts on top of graph
    for (i in 1 : 6) {
      j <- 13 + 16 * (i - 1)
      k <- 1 + 16 * (i - 1)
      text(bp[j], ymax * 1.15, labels = sum(catalog[k : (16 * i), ]), 
           xpd = NA, cex = 0.8)
    }
  }
  
  if (type == "signature") {
    # Calculate mutation signatures of the input catalog
    sig <- catalog / sum(catalog)
    
    # Get ylim
    ymax <- max(sig)
    
    # Barplot
    bp <- barplot(sig[, 1], xaxt = "n", yaxt = 'n', xaxs = "i", xlim = c(-1, 230),
                  lwd = 3, space = 1.35, border = NA, 
                  col = cols, ylab = "proportion")
  }
  
  # Draw the x axis
  Axis(side = 1, at = c(bp[seq(1, 93, 4)] - 0.5, bp[96] + 0.5), 
       labels = FALSE, lwd.tick = 0, lwd = 0.5)
  
  # Draw the y axis
  y.axis.values <- c(0, ymax)
  if (type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  } else {
    y.axis.labels <- y.axis.values
  }
  Axis(side = 2, at = y.axis.values, las = 1, cex.axis = 0.8, labels = FALSE)
  text(-3.5, y.axis.values, labels = y.axis.labels, cex = 0.8,
       las = 1, adj = 1, xpd = NA)
  
  # Draw the ID information on top of graph
  text(bp[2], ymax * 1.08, labels = id, xpd = NA, font = 2, adj = c(0, 0))
  
  # Draw the labels along x axis
  xlabel.idx <- seq(1, 96, by = 4)
  label <- c("A", "C", "G", "T")
  
  # Draw the first line of x axis label
  text(bp[xlabel.idx], -ymax / 7, labels = label, 
       cex = 0.7, adj = 0.5, xpd = NA)
  
  x <- list(bp[xlabel.idx], bp[xlabel.idx + 1], 
            bp[xlabel.idx + 2], bp[xlabel.idx + 3])
  y <- c(-ymax / 3.5, -ymax / 2.8, -ymax / 2.5, -ymax / 2.1)
  # Draw the remaining lines of x axis labels
  for (i in 1 : 4) {
    text(x[[i]], y[i], labels = label[i], cex = 0.7, adj = 0.5, xpd = NA)
  } 
  
  # Draw the text on the left plane
  text(1.5, -ymax / 7, labels = "preceded by 5'", 
       pos = 2, xpd = NA, cex = 0.8)
  text(1.5, -ymax / 3.5, labels = "followed by 3'", 
       pos = 2, xpd = NA, cex = 0.8)
  
  # Draw horizontal lines and names of major mutation class on top of graph
  x.left <- bp[seq(1, 81, 16)] 
  x.right <- bp[seq(16, 96, 16)] 
  rect(xleft = x.left, ymax * 1.28, xright = x.right, ymax * 1.3, 
       col = class.col, border = NA, xpd = NA, adj = 0.5)
  text((x.left + x.right)/2, ymax * 1.38, labels = maj.class.names, xpd = NA)
}

Cat96ToPdf <- 
  function(catalog, name, id, type = "density", abundance = NULL) {
    # Plot the SNS 96 mutation catalog of different samples to a PDF file.
    # 
    # Args:
    #   catalog:   A matrix whose rownames indicate the 96 SNS mutation types 
    #              while its columns contain the counts of each mutation type
    #              from different samples.
    #   name:      The name of the PDF file to be produced.
    #   id:        A vector containing the ID information of different samples.
    #   type:      A vector of values indicating the type of plot for each sample.
    #              If type = "density", the graph will plot the rates of 
    #              mutations per million trinucleotides for each mutation type.
    #              If type = "counts", the graph will plot the occurrences of   
    #              the 96 mutation types in the sample.
    #              If type = "signature", the graph will plot mutation 
    #              signatures of the sample.
    #              The default value for type is "density".
    #   abundance: A matrix containing trinucleotide abundance information.
    #              To be used only when type = "density".
    
    # Setting the width and length for A4 size plotting
    cairo_pdf(name, width = 8.2677, height = 11.6929, onefile = TRUE) 
    
    n <- ncol(catalog)
    par(mfrow = c(8, 1), mar = c(4, 5.5, 2, 1), oma = c(1, 1, 2, 1))
    
    # Do recycling of the function parameters if a vector 
    # with length more than one is not specified by the user.
    if (n > 1 && length(type) == 1) {
      type <- rep(type, n)
    }
    
    for (i in 1 : n) {
      PlotCat96(catalog[, i, drop = FALSE], 
                id = id[i],
                type = type[i],
                abundance = abundance)
    }
    dev.off()
  }

PlotCat192 <- function(catalog, id, type = "counts", 
                       cex = 0.8, abundance = NULL) {
  # Plot the SNS 192 mutation catalog of one sample.
  # 
  # Args:
  #   catalog:   A matrix whose rownames indicate the 192 SNS mutation types
  #              while its column contains the counts of each mutation type.
  #   id:        The ID information of the sample which has mutations.
  #   type:      A value indicating the type of the graph. 
  #              If type = "counts", the graph will plot the occurrences of   
  #              the 192 mutation types in the sample.
  #              If type = "signature", the graph will plot mutation 
  #              signatures of the sample.
  #              If type = "density", the graph will plot the rates of 
  #              mutations per million trinucleotides for each mutation type.
  #              The default value for type is "counts".
  #   cex:       A numerical value giving the amount by which mutation class 
  #              labels on top of graph, y axis labels and sample name 
  #              should be magnified relative to the default.
  #   abundance: A matrix containing trinucleotide abundance and strand 
  #              information, to be used only when type = "density".
  
  stopifnot(dim(catalog) == c(192, 1))
  
  class.col  <- c("#03bcee",
                  "#010101",
                  "#e32926",
                  "#999999",
                  "#a1ce63",
                  "#ebc6c4")
  
  bg.class.col  <- c("#DCF8FF",
                     "#E9E9E9",
                     "#FFC7C7",
                     "#F7F7F7",
                     "#E5F9DF",
                     "#F9E7E7")
  
  strand.col <- c("#394398", 
                  "#e83020")
  
  # Sort data in plotting order
  counts <- catalog[.to.reorder.192.for.plotting, ]
  
  num.classes <- length(counts)
  maj.class.names = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  cols <- rep(strand.col, num.classes / 2)
  
  if (type == "counts") {
    # Get ylim
    ymax <- max(counts) * 1.3
    
    # Barplot: side by side
    mat <- matrix(counts, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), 
                  axes = FALSE, ann = FALSE, lwd = 3, xaxs = "i",
                  border = NA, col = cols, xpd = NA, ylab = "counts")
  }
  
  if (type == "signature") {
    # Calculate mutation signatures of the input catalog
    sig <- counts / sum(counts)
    
    # Get ylim
    ymax <- ifelse(max(sig) * 1.3 > 1, 1, max(sig) * 1.3)
    
    # Barplot: side by side
    mat <- matrix(sig, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), 
                  axes = FALSE, ann = FALSE, lwd = 3, xaxs = "i",
                  border = NA, col = cols, xpd = NA, ylab = "proportion")
  }
  
  if (type == "density") {
    stop("not implemented")
  }
  
  # Draw lines above each class:
  x.left <- bp[seq(0, 160, 32) + 1] - 0.5
  x.right <- bp[seq(32, 192, 32)] + 0.5
  rect(xleft = x.left, ymax * 1.01, xright = x.right, ymax * 1.03, 
       col = class.col, border = NA, xpd = NA)
  
  # Draw mutation class labels at the top of the figure:
  text((x.left + x.right)/2, ymax * 1.07, labels = maj.class.names, 
       cex = cex, xpd = NA)
  
  # Draw background color
  rect(xleft = x.left - 0.5, 0, xright = x.right + 0.5, ymax, 
       col = bg.class.col, border = 'grey90', lwd = 1.5)
  
  # Plot again
  barplot(mat, beside = TRUE, ylim = c(0, ymax), 
          axes = FALSE, ann = FALSE, lwd = 3, 
          border = NA, col = cols, xpd = NA, add = TRUE)
  
  # Draw grid lines
  segments(bp[1] - 1, seq(0, ymax, ymax/4), bp[num.classes] + 1, 
           seq(0, ymax, ymax/4), col = 'grey35', lwd = 0.25)
  
  # Draw y axis and write mutation counts on top of graph(if applicable)
  y.axis.values <- seq(0, ymax, ymax/4)
  if (type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  } else {
    y.axis.labels <- round(y.axis.values, 0)
    
    # Write the mutation counts on top of graph
    for (i in 1 : 6) {
      j <- 23 + 32 * (i - 1)
      k <- 1 + 32 * (i - 1)
      text(bp[j], ymax * 7 / 8, labels = sum(counts[k : (32 * i)]), 
           xpd = NA, cex = 0.8)
    }
  }
  text(-0.5, y.axis.values, labels = y.axis.labels, 
       las = 1, adj = 1, xpd = NA, cex = cex)
  
  # Draw the x axis labels
  context.pos <- (bp[seq(1, 191, 2)] + bp[seq(2, 192, 2)]) / 2
  xlabel.1 <- c("A", "C", "G", "T")
  xlabel.2 <- rep(c("A", "C", "G", "T"), each = 4)
  text(context.pos, -ymax / 100, labels = rep(xlabel.1, 24), cex = 0.5, 
       srt = 90, adj = 1, xpd = NA)
  text(context.pos, -ymax / 18, labels = rep(c("C", "T"), each = 48), 
       cex = 0.5, srt = 90, adj = 1, xpd = NA)
  text(context.pos, -ymax / 10, labels = rep(xlabel.2, 6), 
       cex = 0.5, srt = 90, adj = 1, xpd = NA)
  
  # Write the name of the sample
  text(1.5, ymax * 7 / 8, labels = id, adj = 0, cex = cex, font = 2)
}

Cat192ToPdf <- function(catalog, name, id, type = "counts",
                        cex = 0.8, abundance = NULL) {
  # Plot the SNS 192 mutation catalog of different samples to a PDF file.
  # 
  # Args:
  #   catalog:   A matrix whose rownames indicate the 192 SNS mutation types
  #              while its columns contain the counts of each mutation type
  #              from different samples.
  #   name:      The name of the PDF file to be produced.
  #   id:        The ID information of the sample which has mutations.
  #   type:      A vector of values indicating the type of graph for each sample.
  #              If type = "counts", the graph will plot the occurrences of   
  #              the 192 mutation types in the sample.
  #              If type = "signature", the graph will plot mutation 
  #              signatures of the sample.
  #              If type = "density", the graph will plot the rates of 
  #              mutations per million trinucleotides for each mutation type.
  #              The default value for type is "counts".
  #   cex:       A numerical value giving the amount by which mutation class 
  #              labels on top of graph, y axis labels and sample name 
  #              should be magnified relative to the default.
  #   abundance: A matrix containing trinucleotide abundance and strand 
  #              information, to be used only when type = "density".
  
  # Setting the width and length for A4 size plotting
  cairo_pdf(name, width = 8.2677, height = 11.6929, onefile = TRUE)
  
  n <- ncol(catalog)
  par(mfrow = c(8, 1), mar = c(2, 4, 2, 2), oma = c(3, 2, 1, 1))
  
  # Do recycling of the function parameters if a vector 
  # with length more than one is not specified by the user.
  if (n > 1 && length(type) == 1) {
    type <- rep(type, n)
  }
  
  for (i in 1 : n) {
    PlotCat192(catalog[, i, drop = FALSE], 
               id = colnames(catalog)[i], type = type[i], 
               cex = cex, abundance = abundance)
  }
  dev.off()
}

PlotCat192Strand <- function(catalog, id, type = "counts",
                             cex = 1, abundance = NULL) {
  # Plot the transcription strand bias graph of 6 SNS mutation types
  # ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") in one sample.
  # 
  # Args:
  #   catalog:   A matrix whose rownames indicate the 192 SNS mutation types
  #              while its column contains the counts of each mutation type.
  #   id:        The ID information of the sample which has mutations.
  #   type:      A value indicating the type of the graph. 
  #              If type = "counts", the graph will plot the occurrences of   
  #              the 6 SNS mutation types in the sample.
  #              If type = "signature", the graph will plot mutation 
  #              signatures of the 6 SNS mutation types in the sample.
  #              If type = "density", the graph will plot the rates of 
  #              mutations per million trinucleotides for each of the 
  #              6 SNS mutation types.
  #              The default value for type is "counts".
  #   cex:       A numerical value giving the amount by which mutation class 
  #              labels, y axis labels, sample name and legend 
  #              should be magnified relative to the default.
  #   abundance: A matrix containing trinucleotide abundance and strand 
  #              information, to be used only when type = "density".
  
  stopifnot(dim(catalog) == c(192, 1))
  
  strand.col <- c('#394398', 
                  '#e83020')
  
  # Sort data in plotting order
  counts <- catalog[.to.reorder.192.for.plotting, ]
  
  # Get the counts for each major mutation class
  counts.strand <- integer(12)
  for (i in 1 : 6){
    counts.strand[2 * i - 1] <- 
      sum(counts[seq(32 * (i - 1) + 1, by = 2, length.out = 16)])
    counts.strand[2 * i] <- 
      sum(counts[seq(32 * (i - 1) + 2, by = 2, length.out = 16)])
  }
  
  num.classes <- length(counts.strand)
  maj.class.names <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  cols <- rep(strand.col, num.classes / 2)
  
  if (type == "counts") {
    # Get ylim
    ymax <- max(counts.strand) * 1.3
    
    # Barplot: side by side
    mat <- matrix(counts.strand, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 5.5),
                  width = 0.3, xaxs = "i", yaxs = "i",
                  axes = FALSE, ann = FALSE, ylab = "counts",
                  border = NA, col = cols, xpd = NA)
  }
  
  if (type == "signature") {
    # Calculate mutation signatures of each major mutation class
    sig <- counts.strand / sum(counts.strand)
    
    # Get ylim
    ymax <- ifelse(max(sig) * 1.3 > 1, 1, max(sig) * 1.3)
    
    # Barplot: side by side
    mat <- matrix(sig, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 5.5),
                  width = 0.3, xaxs = "i", yaxs = "i",
                  axes = FALSE, ann = FALSE, ylab = "proportion",
                  border = NA, col = cols, xpd = NA)
  }
  
  if (type == "density") {
    stop("not implemented")
  }
  
  # Draw y axis
  y.axis.values <- seq(0, ymax, length.out = 5)
  if (type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  } else {
    y.axis.labels <- round(y.axis.values, 0)
  }
  Axis(side = 2, at = y.axis.values, las = 1, labels = FALSE)
  text(-0.25, y.axis.values, labels = y.axis.labels, 
       las = 1, adj = 1, xpd = NA, cex = cex)
       
  # Draw x axis
  xlabel <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  text(bp[seq(1, 12, by = 2)] + 0.16, -ymax / 28, 
       labels = xlabel, xpd = NA, cex = cex)
  
  # Add legend
  legend(bp[6], ymax * 0.95, fill = strand.col, border = "white", 
         xpd = NA, bty = "n",
         legend = c("Transcribed", "Untranscribed"), cex = cex)
  
  # Draw the ID information on top of graph
  text(bp[5], ymax * 1.02, labels = id, xpd = NA, 
       font = 2, cex = cex, adj = c(0, 0))
}

Cat192StrandToPdf <- function(catalog, name, id, type = "counts",
                              cex = 1, abundance = NULL) {
  # Plot the transcription strand bias graph of 6 SNS mutation types
  # ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") of different samples to a PDF file.
  #
  # Args:
  #   catalog:   A matrix whose rownames indicate the 192 SNS mutation types
  #              while its columns contain the counts of each mutation type
  #              from different samples.
  #   name:      The name of the PDF file to be produced.
  #   id:        The ID information of the sample which has mutations.
  #   type:      A vector of values indicating the type of graph for each sample.
  #              If type = "counts", the graph will plot the occurrences of   
  #              the 192 mutation types in the sample.
  #              If type = "signature", the graph will plot mutation 
  #              signatures of the sample.
  #              If type = "density", the graph will plot the rates of 
  #              mutations per million trinucleotides for each mutation type.
  #              The default value for type is "counts".  
  #   cex:       A numerical value giving the amount by which mutation class 
  #              labels, y axis labels, sample name and legend 
  #              should be magnified relative to the default.
  #   abundance: A matrix containing trinucleotide abundance and strand 
  #              information, to be used only when type = "density".
  
  # Setting the width and length for A4 size plotting
  cairo_pdf(name, width = 8.2677, height = 11.6929, onefile = TRUE)
  
  n <- ncol(catalog)
  par(mfrow = c(4, 3), mar = c(2, 5, 2, 1), oma = c(2, 2, 2, 2))
  
  # Do recycling of the function parameters if a vector 
  # with length more than one is not specified by the user.
  if (n > 1 && length(type) == 1) {
    type <- rep(type, n)
  }
  
  for (i in 1 : n) {
    PlotCat192Strand(catalog[, i, drop = FALSE], 
               id = colnames(catalog)[i], type = type[i],
               cex = cex, abundance = abundance)
  }
  dev.off()
}

# Setting up functions for PlotCat1536

# Define the bases and their colors in plot
bases <- c("A", "C", "G", "T")
base.cols <- c("forestgreen", "dodgerblue2", "black", "red")

# Define the theme color for plotting
theme.col <- "darkgreen"
colPal <- colorRampPalette(c("white", theme.col))

scale.col <- function(x, x.max) {
  # Do the scaling of colors.
  # 
  # Args:
  #   x:      An integer number.
  #   x.max:  An integer number no less than x.
  #
  # Returns: 
  #   A scaled color type. 
  
  idx <- round(x / x.max * 256)
  col <- colPal(256)[idx]
  return(col)
}

DrawImage <- function(counts, colors) {
  # Plot the color images of the mutation counts, which 
  # has its box, axes and axes' annotations removed.
  # 
  # Args:
  #   counts: The number of mutation counts for pentanucleotides in one 
  #           sample, normalized by pentanucleotide occurrence in the genome.
  #   colors: The colors used in plotting.  
  
  image(1 : 16, 1 : 16, matrix(counts, 16, 16), 
        col = colors, axes = FALSE, ann = FALSE, bty = "n")
  
  # Draw the gridlines
  segments(rep(0.5, 5), c(0.5, 4.5, 8.5, 12.5, 16.5), 
           rep(16.5, 5), c(0.5, 4.5, 8.5, 12.5, 16.5), xpd = T)
  segments(c(0.5, 4.5, 8.5, 12.5, 16.5), rep(0.5, 5), 
           c(0.5, 4.5, 8.5, 12.5, 16.5), rep(16.5, 5), xpd = T)
}

DrawAxisY <- function() {
  # Draw the Y axis of the plot with its annotations and gridlines.
  
  text(0, 16 : 1, rep(bases, each = 4), col = rep(base.cols, each = 4), 
       srt = 90, font = 2, xpd = T)
  text(-1, 16 : 1, rep(bases, 4), col = rep(base.cols, 4), 
       srt = 90, font = 2, xpd = T)
  text(-2.5, 8.5, "Preceding bases", srt = 90, cex = 1.5, xpd = T)
  segments(rep(-1.3, 5), c(0.5, 4.5, 8.5, 12.5, 16.5), rep(0.5, 5), 
           c(0.5, 4.5, 8.5, 12.5, 16.5), xpd = T)
}

DrawAxisX <- function() {
  # Draw the X axis of the plot with its annotations and gridlines.
  
  text(1 : 16, 17, rep(bases, each = 4), col=rep(base.cols, each = 4), 
       font = 2, xpd = T)
  text(1 : 16, 18, rep(bases, 4), col = rep(base.cols, 4), 
       font = 2, xpd = T)
  segments(c(0.5, 4.5, 8.5, 12.5, 16.5), rep(16.5, 5), 
           c(0.5, 4.5, 8.5, 12.5, 16.5), rep(18.3, 5), xpd = T)
}

PlotCat1536 <- 
  function(catalog, id, scale = TRUE, abundance) {
    # Plot the pentanucleotide sequence contexts for one sample, 
    # normalized by pentanucleotide occurrence in the genome.
    # 
    # Args:
    #   catalog:   A matrix whose rownames indicate the 1536 SNS mutation types
    #              while its column contains the counts of each mutation type.
    #              The mutation types are in six-letters like CATTAT, 
    #              first 2-letters CA refers to (-2, -1) position,
    #              third letter T refers to the base which has mutation,
    #              next second 2-letters TA refers to (+1, +2) position,
    #              last letter T refers to the base after mutation.
    #   id:        The id of the sample to be displayed on top of the graph.
    #   scale:     A logical value indicating whether to do color scaling
    #              for all mutation types.
    #   abundance: A matrix containing pentanucleotide abundance information.
    
    stopifnot(dim(catalog) == c(1536, 1))
    mut.type <- rownames(catalog)
    
    # Calculate pentanucleotide sequence, normalized by pentanucleotide 
    # occurrence in the genome
    rates <- catalog
    for (i in 1 : 1536) {
      penta.names <- substr(mut.type[i], 1, 5)
      rates[i] <- 
        catalog[i] * 1000000 / abundance[penta.names, ]
    }
    
    # Sort the rates matrix in plotting order
    rates <- as.data.frame(rates)
    rates$mut.type <- mut.type
    rates$ref2alt <- paste0(substr(mut.type, 3, 3), substr(mut.type, 6, 6))
    rates$minus2bs <- substr(mut.type, 1, 1)
    rates$minus1bs <- substr(mut.type, 2, 2)
    rates$plus1bs <- substr(mut.type, 4, 4)
    rates$plus2bs <- substr(mut.type, 5, 5)
    rates <- arrange(rates, ref2alt, desc(minus1bs), desc(minus2bs), 
                     plus1bs, plus2bs)
    plot.order <- rates$mut.type
    rates <- as.matrix(rates[, 1])
    rownames(rates) <- plot.order
    
    main.mut.type <- 
      paste(substr(plot.order, 3, 3), substr(plot.order, 6, 6), sep = '>')
    main.types <- unique(main.mut.type)
    n.types <- length(main.types)  # max 6 types
    
    # Plot one sample on one page
    par(mfrow = c(2, 3), oma = c(1, 1, 1, 1)) 
    
    # Do the color scaling 
    for (i in 1 : n.types) {
      main.type <- main.types[i]
      sub.rates <- rates[main.mut.type == main.type, 1, drop = FALSE]
      if (scale) {
        max.col <- scale.col(max(sub.rates), max(rates))
      } else {
        max.col <- "darkgreen"
      }
      
      col.ref <- colorRampPalette(c("white", max.col))(256)
      
      # Draw the 6 plots on page one by one
      if (i == 1) {
        par(mar = c(1, 5, 7.5, 1))
        DrawImage(sub.rates, col.ref)
        DrawAxisY()
        DrawAxisX()
        text(8.5, 19, main.type, cex = 1.5, xpd = T)
      }
      
      if (i == 2) {
        par(mar = c(1, 1, 7.5, 1))
        DrawImage(sub.rates, col.ref)
        DrawAxisX()
        text(8.5, 19, main.type, cex = 1.5, xpd = T)
        text(8.5, 20.5, id, cex = 1.5, xpd = T)
      }
      
      if (i == 3) {
        par(mar = c(1, 1, 7.5, 3))
        DrawImage(sub.rates, col.ref)
        DrawAxisX()
        text(8.5, 19, main.type, cex = 1.5, xpd = T)
        text(17.5, 17, '1bp 3\'', xpd = T, cex = 1)
        text(17.5, 18, '2bp 3\'', xpd = T, cex = 1)
      }
      
      if (i == 4) {
        par(mar = c(2.5, 5, 2, 1))
        DrawImage(sub.rates, col.ref)
        DrawAxisY()
        text(8.5, 17, main.type, cex = 1.5, xpd = T)
        text(-1, -0.7, '1bp 5\'', xpd = T, srt = 45, adj = 0, cex = 1)
        text(-2, -0.7, '2bp 5\'', xpd = T, srt = 45, adj = 0, cex = 1)
      }
      
      if (i == 5) {
        par(mar = c(2.5, 1, 2, 1))
        DrawImage(sub.rates, col.ref)
        text(8.5, 17, main.type, cex = 1.5, xpd = T)
      }
      
      if (i == 6) {
        par(mar = c(2.5, 1, 2, 3))
        DrawImage(sub.rates, col.ref)
        text(8.5, 17, main.type, cex = 1.5, xpd = T)
      }
    }
  }

#' Plot the 1536 mutation catalog of >= 1 samples to a PDF file.
#' 
#' @param catalog   A matrix whose rownames indicate the 1536 SNS mutation types
#'              while its columns contain the counts of each mutation type
#'              from different samples.
#'              The mutation types are in six-letters like CATTAT, 
#'              first 2-letters CA refers to (-2, -1) position,
#'              third letter T refers to the base which has mutation,
#'              next second 2-letters TA refers to (+1, +2) position,
#'              last letter T refers to the base after mutation.
#' @param name Name of the PDF file to be produced.
#' @param id  A vector containing the identifier of each sample.
#' @param abundance A matrix containing pentanucleotide abundance information.
Cat1536ToPdf <- function(catalog, name, id, abundance) {
  
  cairo_pdf(name, width = 11.6929, height = 9.2677, onefile = TRUE)
  
  n <- ncol(catalog)
  
  for (i in 1 : n) {
    PlotCat1536(catalog[, i, drop = FALSE], 
                id = colnames(catalog)[i], 
                abundance = abundance)
  }
  dev.off()
}

###############################################################################
# Plotting functions for DNS 78, 144 catalog start here
###############################################################################

PlotCatDNS78 <- function(catalog, id, type = "density", abundance = NULL) {
  # Plot the DNS 78 mutation catalog of one sample.
  # 
  # Args:
  #   catalog:   A matrix whose rownames indicate the 78 DNS mutation types 
  #              while its columns contain the counts of each mutation type
  #              from different samples.
  #   id:        The ID information of the sample which has mutations.
  #   type:      A value indicating the type of the graph. 
  #              If type = "density", the graph will plot the rates of 
  #              mutations per million nucleotides for each mutation type.
  #              If type = "counts", the graph will plot the occurrences of   
  #              the 78 mutation types in the sample.
  #              If type = "signature", the graph will plot mutation 
  #              signatures of the sample.
  #              The default value for type is "density".
  #   xlabel:    A logical value indicating whether to draw the x axis labels.
  #   upper:     A logical value indicating whether to draw horizontal lines 
  #              and names of major mutation class on top of graph.
  #   abundance: A matrix containing dinucleotide abundance information, 
  #              to be used only when type = "density".
  
  stopifnot(dim(catalog) == c(78, 1))
  stopifnot(rownames(catalog) == .catalog.row.order.DNS.78)
  
  dinuc.class.col <- brewer.pal(10, "Paired")
  cols <- rep(dinuc.class.col, c(9, 6, 9, 6, 9, 6, 6, 9, 9, 9))
  num.classes <- length(catalog)
  maj.class.names <- 
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  
  if (type == "density") {
    # Calculate rate of mutations per million nucleotides for the catalog
    rate <- double(78)
    for (i in 1 : 78) {
      rate[i] <- 
        catalog[i] * 1000000 / abundance[substr(rownames(catalog)[i], 1, 2), ]
    }
    
    # Get ylim
    ymax <- ifelse(max(rate) * 1.3 > 1, 1, max(rate) * 1.3)
    
    # Barplot
    bp <- barplot(rate, ylim = c(0, ymax), xaxt = "n", yaxt = "n", xaxs = "i",  
                  lwd = 3, space = 1.35, border = NA, col = cols, 
                  xpd = NA, ylab = "mut/million")
    
    # Write the mutation counts on top of graph
    for (i in 1 : 10) {
      j <- c(7, 13, 20, 28, 35, 43, 49, 56, 65, 74)
      name <- substr(rownames(catalog), 1, 2)
      text(bp[j[i]], ymax * 0.9, xpd = NA, cex = 0.8,
           labels = sum(catalog[name == maj.class.names[i], ]))
    }
  }
  
  if (type == "counts") {
    # Get ylim
    ymax <- max(catalog[, 1]) * 1.3
    
    # Barplot
    bp <- barplot(catalog[, 1], xaxt = "n", yaxt = "n", ylim = c(0, ymax),
                  xaxs = "i", lwd = 3, space = 1.35, border = NA, xaxs = "i",
                  col = cols, ylab = "counts")
    
    # Write the mutation counts on top of graph
    for (i in 1 : 10) {
      j <- c(7, 13, 20, 28, 35, 43, 49, 56, 65, 74)
      name <- substr(rownames(catalog), 1, 2)
      text(bp[j[i]], ymax * 0.9, xpd = NA, cex = 0.8,
           labels = sum(catalog[name == maj.class.names[i], ]))
    }
  }
    
  if (type == "signature") {
    # Calculate mutation signatures of the input catalog
    sig <- catalog / sum(catalog)
    
    # Get ylim
    ymax <- ifelse(max(sig) * 1.3 > 1, 1, max(sig) * 1.3)
    
    # Barplot
    bp <- barplot(sig[, 1], xaxt = "n", yaxt = "n", ylim = c(0, ymax),
                  lwd = 3, space = 1.35, border = NA, xaxs = "i",
                  col = cols, ylab = "proportion")
  }
  
  # Draw box and grid lines
  rect(xleft = bp[1] - 1, 0, xright = bp[num.classes] + 1, ymax, 
       border = "grey60", lwd = 0.5, xpd = NA)
  segments(bp[1] - 1, seq(0, ymax, ymax / 4), bp[num.classes] + 1, 
           seq(0, ymax, ymax / 4), col = "grey60", lwd = 0.5, xpd = NA)

  # Draw lines above each class:
  x.left <- bp[c(0, 9, 15, 24, 30, 39, 45, 51, 60, 69) + 1] - 0.5
  x.right <- bp[c(9, 15, 24, 30, 39, 45, 51, 60, 69, 78)] + 0.5
  rect(xleft = x.left, ymax * 1.01, xright = x.right, ymax * 1.08, 
       col = dinuc.class.col, border = NA, xpd = NA)
  
  # Draw mutation class labels at the top of the graph
  text((x.left + x.right) / 2, ymax * 1.123, 
       labels = paste(maj.class.names, "NN", sep = ">"), cex = 0.7, xpd = NA)
  
  # Draw the ID information on top of graph
  text(1.5, ymax * 7 / 8, labels = id, adj = 0, cex = 0.8, font = 2)
  
  # Draw y axis
  y.axis.values <- seq(0, ymax, ymax / 4)
  if (type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  } else {
    y.axis.labels <- round(y.axis.values, 0)
  }
  text(0.35, y.axis.values, labels = y.axis.labels, 
       las = 1, adj = 1, xpd = NA, cex = 0.75)
  
  # Draw x axis labels
  text(bp, -ymax / 80, labels = substr(rownames(catalog), 4, 4), 
       cex = 0.5, srt = 90, adj = 1, xpd = NA)
  text(bp, -ymax / 15, labels = substr(rownames(catalog), 3, 3), 
       cex = 0.5, srt = 90, adj = 1, xpd = NA)
}

CatDNS78ToPdf <- 
  function(catalog, name, id, type = "density", abundance = NULL) {
    # Plot the DNS 78 mutation catalog of different samples to a PDF file.
    # 
    # Args:
    #   catalog:   A matrix whose rownames indicate the 78 DNS mutation types 
    #              while its columns contain the counts of each mutation type
    #              from different samples.
    #   name:      The name of the PDF file to be produced.
    #   id:        A vector containing the ID information of different samples.
    #   type:      A vector of values indicating the type of plot for each sample.
    #              If type = "density", the graph will plot the rates of 
    #              mutations per million nucleotides for each mutation type.
    #              If type = "counts", the graph will plot the occurrences of   
    #              the 78 mutation types in the sample.
    #              If type = "signature", the graph will plot mutation 
    #              signatures of the sample.
    #              The default value for type is "density".
    #   abundance: A matrix containing dinucleotide abundance information,  
    #              to be used only when type = "density".
    
    # Setting the width and length for A4 size plotting
    cairo_pdf(name, width = 8.2677, height = 11.6929, onefile = TRUE) 
    
    n <- ncol(catalog)
    par(mfrow = c(8, 1), mar = c(2, 4, 2, 2), oma = c(3, 3, 2, 2))
    
    # Do recycling of the function parameters if a vector 
    # with length more than one is not specified by the user.
    if (n > 1 && length(type) == 1) {
      type <- rep(type, n)
    }
    
    for (i in 1 : n) {
      PlotCatDNS78(catalog[, i, drop = FALSE], 
                id = id[i],
                type = type[i],
                abundance = abundance)
    }
    dev.off()
  }

PlotCatDNS144 <- function(catalog, id, type = "counts",
                          cex = 1, abundance = NULL) {
  # Plot the transcription strand bias graph of 10 major DNS mutation types
  # ("AC>NN", "AT>NN", "CC>NN", "CG>NN", "CT>NN", "GC>NN", "TA>NN", 
  # "TC>NN", "TG>NN", "TT>NN") in one sample.
  # 
  # Args:
  #   catalog:   A matrix whose rownames indicate the 144 DNS mutation types
  #              while its column contains the counts of each mutation type.
  #   id:        The ID information of the sample which has mutations.
  #   type:      A value indicating the type of the graph. 
  #              If type = "counts", the graph will plot the occurrences of   
  #              the 10 major DNS mutation types in the sample.
  #              If type = "signature", the graph will plot mutation 
  #              signatures of the 10 major DNS mutation types in the sample.
  #              If type = "density", the graph will plot the rates of 
  #              mutations per million dinucleotides for each of the 
  #              10 major DNS mutation types.
  #              The default value for type is "counts".
  #   cex:       A numerical value giving the amount by which mutation class 
  #              labels, y axis labels, sample name and legend 
  #              should be magnified relative to the default.
  #   abundance: A matrix containing dinucleotide abundance and strand 
  #              information, to be used only when type = "density".
  
  stopifnot(dim(catalog) == c(144, 1))
  
  strand.col <- c('#394398', 
                  '#e83020')
  
  # Sort data in plotting order
  counts <- catalog[.to.reorder.144.for.plotting, ]
  
  # Get the counts for each major mutation class
  counts.strand <- integer(20)
  for (i in 1 : 10){
    idx <- c(0, 18, 24, 42, 48, 66, 72, 78, 96, 114, 132)
    counts.strand[2 * i - 1] <- 
      sum(counts[seq(idx[i] + 1, idx[i + 1], by = 2)])
    counts.strand[2 * i] <- 
      sum(counts[seq(idx[i] + 2, idx[i + 1], by = 2)])
  }
  
  num.classes <- length(counts.strand)
  maj.class.names <- 
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  cols <- rep(strand.col, num.classes / 2)
  
  if (type == "counts") {
    # Get ylim
    ymax <- max(counts.strand) * 1.3
    
    # Barplot: side by side
    mat <- matrix(counts.strand, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 9),
                  width = 0.3, xaxs = "i", yaxs = "i",
                  axes = FALSE, ann = FALSE, ylab = "counts",
                  border = NA, col = cols, xpd = NA)
  }
  
  if (type == "signature") {
    # Calculate mutation signatures of each major mutation class
    sig <- counts.strand / sum(counts.strand)
    
    # Get ylim
    ymax <- ifelse(max(sig) * 1.3 > 1, 1, max(sig) * 1.3)
    
    # Barplot: side by side
    mat <- matrix(sig, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 9),
                  width = 0.3, xaxs = "i", yaxs = "i",
                  axes = FALSE, ann = FALSE, ylab = "proportion",
                  border = NA, col = cols, xpd = NA)
  }
  
  if (type == "density") {
    stop("not implemented")
  }
  
  # Draw y axis
  y.axis.values <- seq(0, ymax, length.out = 5)
  if (type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  } else {
    y.axis.labels <- round(y.axis.values, 0)
  }
  Axis(side = 2, at = y.axis.values, las = 1, labels = FALSE)
  text(-0.35, y.axis.values, labels = y.axis.labels, 
       las = 1, adj = 1, xpd = NA, cex = cex)
  
  # Draw x axis
  xlabel <- c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  text(bp[seq(1, 20, by = 2)] + 0.16, -ymax / 28, 
       labels = xlabel, xpd = NA, cex = cex)
  
  # Add legend
  legend(bp[10], ymax * 0.95, fill = strand.col, border = "white", 
         xpd = NA, bty = "n",
         legend = c("Transcribed", "Untranscribed"), cex = cex)
  
  # Draw the ID information on top of graph
  text(bp[8], ymax, labels = id, xpd = NA, 
       font = 2, cex = cex, adj = c(0, 0))
}

CatDNS144ToPdf <- function(catalog, name, id, type = "counts",
                           cex = 1, abundance = NULL) {
  # Plot the transcription strand bias graph of 10 major DNS mutation types
  # ("AC>NN", "AT>NN", "CC>NN", "CG>NN", "CT>NN", "GC>NN", "TA>NN", 
  # "TC>NN", "TG>NN", "TT>NN") of different samples to a PDF file.
  #
  # Args:
  #   catalog:   A matrix whose rownames indicate the 144 DNS mutation types
  #              while its columns contain the counts of each mutation type
  #              from different samples.
  #   name:      The name of the PDF file to be produced.
  #   id:        The ID information of the sample which has mutations.
  #   type:      A vector of values indicating the type of graph for each sample.
  #              If type = "counts", the graph will plot the occurrences of   
  #              the 10 major DNS mutation types in the sample.
  #              If type = "signature", the graph will plot mutation 
  #              signatures of the 10 major DNS mutation types in the sample.
  #              If type = "density", the graph will plot the rates of 
  #              mutations per million dinucleotides for each of the 
  #              10 major DNS mutation types.
  #              The default value for type is "counts".  
  #   cex:       A numerical value giving the amount by which mutation class 
  #              labels, y axis labels, sample name and legend 
  #              should be magnified relative to the default.
  #   abundance: A matrix containing dinucleotide abundance and strand 
  #              information, to be used only when type = "density".
  
  # Setting the width and length for A4 size plotting
  cairo_pdf(name, width = 8.2677, height = 11.6929, onefile = TRUE)
  
  n <- ncol(catalog)
  par(mfrow = c(4, 3), mar = c(2, 5, 2, 1), oma = c(2, 2, 2, 2))
  
  # Do recycling of the function parameters if a vector 
  # with length more than one is not specified by the user.
  if (n > 1 && length(type) == 1) {
    type <- rep(type, n)
  }
  
  for (i in 1 : n) {
    PlotCatDNS144(catalog[, i, drop = FALSE], 
                  id = colnames(catalog)[i], type = type[i],
                  cex = cex, abundance = abundance)
  }
  dev.off()
}

###############################################################################
# Plotting functions for insertion and deletion catalog start here
###############################################################################

PlotCatID <- function(catalog, id, type = "counts"){
  # Plot the insertion and deletion catalog of one sample.
  # 
  # Args:
  #   catalog: A matrix whose rownames indicate the insertion and 
  #            deletion mutation types while its column contains
  #            the counts of each mutation type.
  #   id:      The ID information of the sample which has mutations.
  #   type:    A value indicating the type of the graph. 
  #            If type = "counts", the graph will plot the occurrences of   
  #            the insertion and deletion mutation types in the sample.
  #            If type = "signature", the graph will plot mutation 
  #            signatures of the sample.
  #            The default value for type is "counts".
  
  stopifnot(dim(catalog) == c(83, 1))
  
  indel.class.col <- c("#fdbe6f",
                       "#ff8001",
                       "#b0dd8b",
                       "#36a12e",
                       "#fdcab5",
                       "#fc8a6a",
                       "#f14432",
                       "#bc141a",
                       "#d0e1f2",
                       "#94c4df",
                       "#4a98c9",
                       "#1764ab",
                       "#e2e2ef",
                       "#b6b6d8",
                       "#8683bd",
                       "#61409b")
  
  num.classes <- length(catalog)
  cols <- rep(indel.class.col, 
              c(6, 6, 6, 6,
                6, 6, 6, 6,
                6, 6, 6, 6,
                1, 2, 3, 5))
  
  if (type == "counts") {
    # Get ylim
    ymax <- max(catalog) * 1.3
    
    # Barplot
    bp <- barplot(catalog[, 1], ylim = c(0, ymax), axes = FALSE, xaxt = "n",
                  lwd = 3, space = 1.35, border = NA, col = cols, xpd = NA,
                  xaxs = "i", ylab = "counts")
  }
  
  if (type == "signature") {
    # Calculate mutation signatures of the input catalog
    sig <- catalog / sum(catalog)
    
    # Get ylim
    ymax <- ifelse(max(sig) * 1.3 > 1, 1, max(sig) * 1.3)
    
    # Barplot
    bp <- barplot(sig[, 1], ylim = c(0, ymax), axes = FALSE, xaxt = "n",
                  lwd = 3, space = 1.35, border = NA, col = cols, xpd = NA,
                  xaxs = "i", ylab = "proportion")
  }
  
  # Draw box and grid lines
  rect(xleft = bp[1] - 1.5, 0, xright = bp[num.classes] + 1, ymax, col = NA,
       border = "grey60", lwd = 0.5, xpd = NA)
  segments(bp[1] - 1.5, seq(0, ymax, ymax / 4), bp[num.classes] + 1, 
           seq(0, ymax, ymax / 4), col = "grey60", lwd = 0.5, xpd = NA)
  
  # Draw mutation class labels and lines above each class
  maj.class.names <- c("1bp deletion", "1bp insertion", 
                       ">1bp deletions at repeats\n(Deletion length)",
                       ">1bp insertions at repeats\n(Insertion length)", 
                       "Deletions with microhomology\n(Deletion length)")
  x.left <- bp[c(seq(0, 66, 6), 72, 73, 75, 78) + 1] - 0.5
  x.right <- bp[c(seq(6, 72, 6), 73, 75, 78, 83)] + 0.5
  class.pos <- c((x.left[seq(1, 4, 2)] + x.right[seq(2, 5, 2)]) / 2,
                 (x.left[c(5, 9)] + x.right[c(8, 12)] - 12) / 2, 
                 (x.left[13] + x.right[length(x.left)]) / 2)
  category.lab <- c(rep(c("C", "T"), 2), rep(c("2", "3", "4", "5+"), 3)) 
  category.col <- c(rep(c("black", "white"), 2), 
                    rep(c("black", "black", "black", "white"), 3)) 
  
  # Draw lines above each class
  rect(xleft = x.left, ymax * 1.01, xright = x.right, ymax * 1.09, 
       col = indel.class.col, border = NA, xpd = NA)
  text((x.left + x.right) / 2, ymax * 1.05, labels = category.lab, 
       cex = 0.55, col = category.col, xpd = NA)
  
  # Draw mutation class labels at the top of the figure
  text(class.pos, ymax * 1.23, labels = maj.class.names, cex = 0.75, xpd = NA)
  
  # Draw the ID information of the sample
  text(1.5, ymax * 7 / 8, labels = id, adj = 0, cex = 0.85, font = 2)
  
  # Draw y axis
  y.axis.values <- seq(0, ymax, ymax / 4)
  if (type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  } else {
    y.axis.labels <- round(y.axis.values, 0)
  }
  text(0, y.axis.values, labels = y.axis.labels, 
       las = 1, adj = 1, xpd = NA, cex = 0.75)
  
  # Draw x axis labels
  mut.type <- c(rep(c("1", "2", "3", "4", "5", "6+"), 2), 
                rep(c("0", "1", "2", "3", "4", "5+"), 2), 
                rep(c("1", "2", "3", "4", "5", "6+"), 4), 
                rep(c("0", "1", "2", "3", "4", "5+"), 4), 
                "1", "1", "2", "1", "2", "3", "1", "2", "3", "4", "5+")
  bottom.pos <- c((x.left[1] + x.right[2]) / 2, (x.left[3] + x.right[4]) / 2,
                  class.pos[3 : length(class.pos)])
  bottom.lab <- c("Homopolymer length", "Homopolymer length", 
                  "Number of repeat units", "Number of repeat units",
                  "Microhomology length")
  rect(xleft = x.left, -ymax * 0.09, xright = x.right, -ymax * 0.01, 
       col = indel.class.col, border = NA, xpd = NA)
  text(bp, -ymax * 0.15, labels = mut.type, cex = 0.65, xpd = NA)
  text(bottom.pos, -ymax * 0.27, labels = bottom.lab, cex = 0.75, xpd = NA)
}

CatIDToPdf <- 
  function(catalog, name, id, type = "counts") {
    # Plot the insertion and deletion catalog of different samples to a PDF file.
    # 
    # Args:
    #   catalog: A matrix whose rownames indicate the insertion and 
    #            deletion mutation types while its column contains
    #            the counts of each mutation type from different samples.
    #   name:    The name of the PDF file to be produced.
    #   id:      A vector containing the ID information of different samples.
    #   type:    A vector of values indicating the type of plot for each sample.
    #            If type = "counts", the graph will plot the occurrences of   
    #            the insertion and deletion mutation types in the sample.
    #            If type = "signature", the graph will plot mutation 
    #            signatures of the sample.
    #            The default value for type is "counts".
    
    # Setting the width and length for A4 size plotting
    cairo_pdf(name, width = 8.2677, height = 11.6929, onefile = TRUE) 
    
    n <- ncol(catalog)
    par(mfrow = c(8, 1), mar = c(3, 4, 2.5, 2), oma = c(3, 3, 2, 2))
    
    # Do recycling of the function parameters if a vector 
    # with length more than one is not specified by the user.
    if (n > 1 && length(type) == 1) {
      type <- rep(type, n)
    }
    
    for (i in 1 : n) {
      PlotCatID(catalog[, i, drop = FALSE], 
                id = id[i],
                type = type[i])
    }
    dev.off()
  }
###############################################################################
# Test functions for plotting start here
###############################################################################

TestPlotCat96 <- function() {
  # This function is to test plotting the input SNS 96 mutation catalog. 
  
  par(mar = c(6.5, 5.1, 4.6, 0))
  catalog <- ReadCat96("data/regress.cat.96.csv")
  PlotCat96(catalog[, 1, drop = FALSE], id = "test",  
            type = "density", abundance = .abundance.3bp)
}

TestCat96ToPdf <- function() {
  # This function is to test plotting SNS 96 mutation catalog of 
  # different samples to a PDF file.
  
  catalog <- ReadCat96("data/regress.cat.96.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat2 <- catalog[, 2, drop = FALSE]
  cat3 <- catalog[, 3, drop = FALSE]
  cat4 <- catalog[, 4, drop = FALSE]
  cat <- cbind(cat1, cat1, cat1, cat2, cat2, cat2, cat3, cat3, cat3,
               cat4, cat4, cat4)
  type <- c("counts", "signature", "density")
  Cat96ToPdf(cat, "PlotCat96.test.pdf", id = colnames(cat),
             type = rep(type, 4), xlabel = TRUE, upper = TRUE,
             abundance = .abundance.3bp)
}

TestPlotCat192 <- function() {
  # This function is to test plotting the input SNS 192 mutation catalog. 
  
  par(mar = c(3, 4, 2, 1))
  catalog <- ReadCat192("data/regress.cat.192.csv")
  PlotCat192(catalog[, 1, drop = FALSE], "test")
}

TestPlotCat192Strand <- function() {
  # This function is to test plotting the transcription strand bias graph of 
  # 6 SNS mutation types("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") in one sample.
  
  par(mar = c(5, 8, 5, 1))
  catalog <- ReadCat192("data/regress.cat.192.csv")
  PlotCat192Strand(catalog[, 1, drop = FALSE], "test")
}

TestCat192ToPdf <- function() {
  # This function is to test plotting SNS 192 mutation catalog of
  # different samples to a PDF file.
  
  catalog <- ReadCat192("data/regress.cat.192.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat2 <- catalog[, 2, drop = FALSE]
  cat3 <- catalog[, 3, drop = FALSE]
  cat4 <- catalog[, 4, drop = FALSE]
  cat <- cbind(cat1, cat1, cat2, cat2, cat3, cat3, cat4, cat4)
  type <- c("counts", "signature")
  Cat192ToPdf(cat, "PlotCat192.test.pdf", id = colnames(cat),
              type = rep(type, 4))
}

TestCat192StrandToPdf <- function() {
  # This function is to test plotting the transcription strand bias graph of 
  # 6 SNS mutation types("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") of 
  # different samples to a PDF file.
  
  catalog <- ReadCat192("data/regress.cat.192.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat2 <- catalog[, 2, drop = FALSE]
  cat3 <- catalog[, 3, drop = FALSE]
  cat4 <- catalog[, 4, drop = FALSE]
  cat <- cbind(cat1, cat1, cat1, cat2, cat2, cat2, cat3, cat3, cat3,
               cat4, cat4, cat4)
  type <- rep(c("counts", "signature"), each = 3)
  Cat192StrandToPdf(cat, "PlotCat192Strand.test.pdf", id = colnames(cat),
                    type = rep(type, 2))
}

TestPlotCat1536 <- function() {
  # This function is to test plotting the input SNS 1536 mutation catalog. 
  
  catalog <- ReadCat1536("data/regress.cat.1536.csv")
  PlotCat1536(catalog[, 1, drop = FALSE], "test", abundance = .abundance.5bp)
}

TestCat1536toPdf <- function() {
  # This function is to test plotting SNS 1536 mutation catalog of 
  # different samples to a PDF file.
  
  catalog <- ReadCat1536("data/regress.cat.1536.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  Cat1536ToPdf(catalog, "PlotCat1536.test.pdf",
               id = colnames(catalog1536), 
               abundance = .abundance.5bp)
}

TestPlotCatDNS78 <- function() {
  # This function is to test plotting the input SNS 96 mutation catalog. 
  
  par(mar = c(2, 4, 3, 2))
  catalog <- ReadCatDNS78("data/regress.cat.dns.78.csv")
  PlotCatDNS78(catalog[, 1, drop = FALSE], id = "test",  
            type = "density", abundance = .abundance.2bp)
}

TestCatDNS78ToPdf <- function() {
  # This function is to test plotting DNS 78 mutation catalog of 
  # different samples to a PDF file.
  
  catalog <- ReadCatDNS78("data/regress.cat.dns.78.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat2 <- catalog[, 2, drop = FALSE]
  cat3 <- catalog[, 3, drop = FALSE]
  cat4 <- catalog[, 4, drop = FALSE]
  cat <- cbind(cat1, cat1, cat1, cat2, cat2, cat2, cat3, cat3, cat3,
               cat4, cat4, cat4)
  type <- c("counts", "signature", "density")
  CatDNS78ToPdf(cat, "PlotCatDNS78.test.pdf", id = colnames(cat),
             type = rep(type, 4),
             abundance = .abundance.2bp)
}

TestPlotCatDNS144 <- function() {
  # This function is to test plotting the transcription strand bias graph of 
  # 10 major DNS mutation types("AC>NN", "AT>NN", "CC>NN", "CG>NN", "CT>NN", "GC>NN", 
  # "TA>NN", "TC>NN", "TG>NN", "TT>NN") in one sample.
  
  par(mar = c(5, 8, 5, 1))
  catalog <- ReadCatDNS144("data/regress.cat.dns.144.csv")
  PlotCatDNS144(catalog[, 1, drop = FALSE], "test")
}

TestCatDNS144ToPdf <- function() {
  # This function is to test plotting the transcription strand bias graph of 
  # 10 major DNS mutation types("AC>NN", "AT>NN", "CC>NN", "CG>NN", "CT>NN", "GC>NN", 
  # "TA>NN", "TC>NN", "TG>NN", "TT>NN") of different samples to a PDF file.
  
  catalog <- ReadCatDNS144("data/regress.cat.dns.144.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat2 <- catalog[, 2, drop = FALSE]
  cat3 <- catalog[, 3, drop = FALSE]
  cat4 <- catalog[, 4, drop = FALSE]
  cat <- cbind(cat1, cat1, cat1, cat2, cat2, cat2, cat3, cat3, cat3,
               cat4, cat4, cat4)
  type <- rep(c("counts", "signature"), each = 3)
  CatDNS144ToPdf(cat, "PlotCatDNS144.test.pdf", id = colnames(cat),
                 type = rep(type, 2))
}

TestPlotCatID <- function() {
  # This function is to test plotting the input insertion and deletion 
  # mutation catalog. 
  
  par(mar = c(6, 4, 6, 3))
  catalog <- ReadCatID("data/BTSG_WGS_PCAWG.indels.csv")
  PlotCatID(catalog[, 1, drop = FALSE], id = "test")
}

TestCatIDToPdf <- function() {
  # This function is to test plotting insertion and deletion 
  # mutation catalog of different samples to a PDF file.
  
  catalog <- ReadCatID("data/BTSG_WGS_PCAWG.indels.csv")
  colnames(catalog) <- paste0("Biliary-AdenoCA", 1 : 35)
  cat1 <- catalog[, 1, drop = FALSE]
  cat2 <- catalog[, 2, drop = FALSE]
  cat3 <- catalog[, 3, drop = FALSE]
  cat4 <- catalog[, 4, drop = FALSE]
  cat <- cbind(cat1, cat1, cat2, cat2, cat3, cat3, cat4, cat4)
  type <- c("counts", "signature")
  CatIDToPdf(cat, "PlotCatID.test.pdf", id = colnames(cat),
             type = rep(type, 4))
}