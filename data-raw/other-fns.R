


# TODO(steve): move this to functionality to ICAMS plotting
# Call Cat96ToPdf, with the argument "id" set to the column names of the catalog
Cat96ToPdf0 <- function(catalog, name, type = "signature") {
  invisible(
    Cat96ToPdf(
      catalog = catalog,
      name = name,
      type = type, id = colnames(catalog)))
}

