


#' Removing missing genotypes/lines/entries/GIDs
#'
#' @description This function removes genotyped individuals from a marker matrix that do not meet
#' the minimum threshold of marker information
#'
#' @param X A marker matrix, usually obtained from the routines in FPGC
#' @param threshold A threshold of missing marker information for genotypes in percentage
#' @param plot.missing Logical. Indicating if a histograms of missing information printed
#'
#' @export

rmMissingGenotypes = function(X = NULL, threshold = 50, plot.missing = FALSE) {
  if (sum(is.na(X)) == 0) {
    stop("There are no missing values in X. Have they been imputed?")
  }

  per.missing = (rowSums(is.na(X)) / ncol(X)) * 100  # Computes % of missing genotype data
  X = X[per.missing < threshold, , drop = FALSE]  # Keep genotypes with less than missing info

  if (plot.missing) {
    hist(per.missing, col = "black", border = "white",
         main = "Histogram of missing genotypes",
         xlab = "Percentage of missing data", ylab = "Number of genotypes")
  }

  print(dim(X))
  return(X)
}
