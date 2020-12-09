#### Functions to estimate heritability ####

#' Broad sense heritability in G + E + GxE models
#'
#' @description h2 is a function to estimate broad sense heritaiblity of multienvironmental trials from the variance components of an lme4 model object.
#' It only considers the two-way interaction between genotypes and environments.
#' @param lm Is an lme4 model object with a GxE term.
#' @param varg Is the genetic variance as expressed in the lme4 model object
#' @param varge is the GxE variance as expressed in the lme4 object
#' @param reps Is the name of the column in the data frame that contains the reps
#' @param envs.col Is the name of the column in the data frame that contains the envrionments, treatments, etc.
#' @export
#' @examples
#' H2 = h2(lm = lm, varg = "g", varge = "g:trt", reps = "reps", envs.col = "sites")

h2 = function(lm=lm, varg="g", varge="g:trt", reps=NULL, envs.col = NULL){
    varcomps = VarCorr(lm)
    varcomps = as.data.frame(print(varcomps, comp=c("Variance"))); varcomps=varcomps[,c(1,4,5)]
    var.g = varcomps$vcov[varcomps$grp==varg]
    var.gxe = varcomps$vcov[varcomps$grp==varge]
    var.res = varcomps$vcov[varcomps$grp=="Residual"]
    col.rep = which(names(lm@frame) == reps)
    n.reps = nlevels(lm@frame[,col.rep])
    env.col = which(names(lm@frame) == envs.col)
    n.env = nlevels(lm@frame[,envs.col])
    var.ph = (varcomps$vcov[varcomps$grp==varg] +
                 (varcomps$vcov[varcomps$grp==varge]/n.env) +
                 varcomps$vcov[varcomps$grp=="Residual"]/(n.reps*n.env))
    h.2=var.g/var.ph
    h2 = as.data.frame(cbind(h.2, var.g, var.gxe, var.res, var.ph, n.reps, n.env))
    print(h2)
}

#' Broad sense heritability in one way models
#'
#' @description  h2.single Isa function to calculate the broad sense heritability from an lme4 model object that does not contain GxE interaction
#' It only considers the two-way interaction between genotypes and environments.
#' @param lm Is an lme4 model object with a GxE term.
#' @param reps Is the name of the column in the data frame that contains the reps
#' @export
#' @examples
#' H2 = h2.single(lm = lm, varg = "g", reps = "reps")

h2single = function (lm=lm, varg="g", reps=NULL) {
    varcomps = VarCorr(lm)
    varcomps = as.data.frame(print(varcomps, comp=c("Variance"))); varcomps=varcomps[,c(1,4,5)]
    var.g = varcomps$vcov[varcomps$grp==varg]
    var.res = varcomps$vcov[varcomps$grp=="Residual"]
    col.rep = which(names(lm@frame) == reps)
    n.reps = nlevels(lm@frame[,col.rep])
    var.ph = (varcomps$vcov[varcomps$grp==varg] +
                 varcomps$vcov[varcomps$grp=="Residual"]/n.reps)
    h.2.single=var.g/var.ph
    h2.single = as.data.frame(cbind(h.2.single, var.g, var.res, var.ph, n.reps))
    print(h2.single)
}
