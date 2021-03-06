##%######################################################%##
#                                                          #
####      FUNCTION TO EXTRACT VARIANCE COMPOENENTS      ####
####              ONLY FROM AN LME4 MODEL               ####
#                                                          #
##%######################################################%##

#' Extract variance components
#'
#' @description  Extraction of variance components from an lme4 object
#' @param lm Is a lme4 model objeject
#' @export


varcomps = function(lm) {
  message("EXTRACTING VARIANCE COMPONENTS")
  varcomps = lme4::VarCorr(lm)
  varcomps = as.data.frame(varcomps, comp=c("Variance")); varcomps=varcomps[,c(1,4,5)]
  varcomps = as.data.frame(t(varcomps))
  colnames(varcomps) = varcomps[1 ,]
  varcomps = varcomps["vcov" ,]
  varcomps[1:ncol(varcomps)] = lapply(X = varcomps[1:ncol(varcomps)], FUN = as.numeric)
  varcomps = round(varcomps, 5)
  varcomps
  }
