#' Process lists of lme4 models
#'
#' @description Function suited to prcess lme4 model objects to calculate the broad sense heritability, and statistical parameter of field experiments,
#' i.e, CV, ASED, LSD, etc.,
#' @param lm Is an lme4 model object. It is a model where the genotypes are included as random effect to extract the G variance
#' @param varg Is the genetic variance as expressed in the lme4 model object
#' @param varge is the GxE variance as expressed in the lme4 object
#' @param reps Is the name of the column in the data frame that contains the reps
#' @param envs.col Is the name of the column in the data frame that contains the envrionments, treatments, etc.
#' @param means.gxe Logical to indicate if the GxE adjusted means will be estimated
#' @param which.means A character to indicate which adusted means will be extracted from the model.
#' @param ASED Logical to indicate if the ASED will be calculated. Default is FALSE, to speed computation since in this implementation ASED requires
#' estimation of all pairwise comparisons of means which can be slow if the data frame is large
#' @param gxe.model Logical to indicate if the lme4 object model is a GxE type
#' @param save.output Logical to indicate if the outuput of adjuste means and trial stats will be saved in the computer
#' @param file.name Character to indicate the name of the file if save.output = TRUE.

#' @export


trialstatslist = function(lm=lm, varg="g", varge=NULL,
                       reps=NULL, envs.col = NULL, means.gxe = F, which.means = NULL, ASED = F,
                       gxe.model=F, save.output = T, file.name= "line-means"){
    message("Loading H2 functions")
    #source("C:/Users/LCRESPO/Documents/CIMMYT/Manuscripts/H2 functions.R")
    message("Loading 'lmerTest' and 'stringr' packages")
    #library(lmerTest)
    #library(stringr)

    if (gxe.model == F) {
        message("EXTRACTING VARIANCE COMPONENTS AND COMPUTING HERITABILITY")
        sts = h2single(lm, varg, reps)
        sts = t(sts)
        colnames(sts)="Estimate"
        Grand.mean = as.data.frame(lme4::fixef(lm)[1])
        rownames(Grand.mean)="Grand mean"
        colnames(Grand.mean)="Estimate"

        varcomps = lme4::VarCorr(lm)
        varcomps = as.data.frame(varcomps, comp=c("Variance")); varcomps=varcomps[,c(1,4,5)]

        # fitting the lm as fixed effect to obtain residual DF
        form = summary(lm)$call$formula
        form = deparse(form)

        form.fixed = stringr::str_replace_all(string = form,
                                     pattern = c("\\(1" = "", "\\|" = "", "\\)" = ""))
        form.fixed = stats::formula(paste(form.fixed, collapse = " "))
        #form.fixed = update.formula(form.fixed, value ~ . )
        mt2 = lm(form.fixed, data = lm@frame)
        an = as.data.frame(stats::anova(mt2))
        dof.res = an[nrow(an), 1]

        # fitting the lm with varg as fixed only
        form2 = gsub(pattern = " ", replacement = "", x = form)
        torep = paste("\\(1.\\|.", varg,"\\)", sep = "")
        form.fixed2 = stringr::str_replace(string = form, pattern = torep, replacement = varg)
        form.fixed2 = stats::formula(paste(form.fixed2, collapse = " "))

        mt3 = (lme4::lmer(form.fixed2, data = lm@frame))

        message("COMPUTING LSMEANS AND PAIRWISE COMPARISONS TO CALCULATE ASED")


        line.means = lmerTest::lsmeansLT(model = mt3, which = varg, pairwise = F)

        if(ASED == T){
        line.means.contrasts = suppressMessages(lmerTest::lsmeansLT(model = mt3,
                                                          which = varg, pairwise = T))
        ASED = as.data.frame(mean(line.means.contrasts$'Std. Error'))
                rownames(ASED)="ASED"
        colnames(ASED)= "Estimate"

        CV <- 100*ASED/abs(Grand.mean)
        rownames(CV)="CV"
        colnames(CV)= "Estimate"
        LSD <- ASED*stats::qt(1-0.05/2, dof.res)
        rownames(LSD) = "LSD"
        colnames(LSD)= "Estimate"
        sts = rbind(sts, Grand.mean, CV, LSD, ASED)
        }



        if (save.output == T){

            message("Saving line lsmeans  and stas to folder in csv format")
            nm = paste(file.name, "_MEANS.csv", sep = "")
            nm.sts = paste(file.name, "_STATS.csv", sep = "")
            #line.means = cbind(file.name, line.means)
            utils::write.csv(line.means, file = nm)
            utils::write.csv(sts, file = nm.sts)
        }

        means.sts.list = list(medias = as.data.frame(line.means), sts = sts)
        return(means.sts.list)
        print(sts)
    }

    else {
        message("EXTRACTING VARIANCE COMPONENTS AND COMPUTING HERITABILITY")
        sts = h2(lm, varg,varge, reps, envs.col)
        sts = t(sts)
        colnames(sts)="Estimate"
        Grand.mean = as.data.frame(lme4::fixef(lm)[1])
        rownames(Grand.mean)="Grand mean"
        colnames(Grand.mean)="Estimate"

        varcomps = lme4::VarCorr(lm)
        varcomps = as.data.frame(varcomps, comp=c("Variance")); varcomps=varcomps[,c(1,4,5)]

        # fitting the lm as fixed effect to obtain residual DF
        form = summary(lm)$call$formula
        form = deparse(form)

        form.fixed = stringr::str_replace_all(string = form,
                                     pattern = c("\\(1" = "", "\\|" = "", "\\)" = ""))
        form.fixed = stats::formula(paste(form.fixed, collapse = " "))
        mt2 = lm(form.fixed, data = lm@frame)
        an = as.data.frame(stats::anova(mt2))
        dof.res = an[nrow(an), 1]

        # fitting the lm with varg as fixed only
        message("FITTING THE FIXED MODEL")

        form2 = gsub(pattern = " ", replacement = "", x = form)
        torep = paste("\\(1.\\|.", varg,"\\)", sep = "")
        form.fixed2 = stringr::str_replace(string = form, pattern = torep, replacement = varg)
        torep2 = paste("\\(1.\\|.",varge,"\\)", sep = "")
        form.fixed3 = stringr::str_replace(string = form.fixed2, pattern = torep2, replacement = varge)
        form.fixed3 = stats::formula(paste(form.fixed3, collapse = " "))
        mt3 = (lme4::lmer(form.fixed3, data = lm@frame))

        message("COMPUTING LSMEANS AND PAIRWISE COMPARISONS TO CALCULATE ASED")


        if(means.gxe == T) {
            line.means = (lmerTest::lsmeansLT(model = mt3, which = which.means, pairwise = F))

        }
        else {
        line.means = suppressMessages(lmerTest::lsmeansLT(model = mt3, which = varg, pairwise = F))
        }

        if(ASED == T) {
        line.means.contrasts = suppressMessages(lmerTest::lsmeansLT(model = mt3,
                                                          which = varg, pairwise = T))
        #SED = as.data.frame(contrast(line.means, method = "pairwise"))
        ASED = as.data.frame(mean(line.means.contrasts$'Std. Error'))
        rownames(ASED)="ASED"
        colnames(ASED)= "Estimate"


        CV <- 100*ASED/abs(Grand.mean)
        rownames(CV)="CV"
        colnames(CV)= "Estimate"
        LSD <- ASED*stats::qt(1-0.05/2, dof.res)
        rownames(LSD) = "LSD"
        colnames(LSD)= "Estimate"
        sts = rbind(sts, Grand.mean, CV, LSD, ASED)
        }

        if (save.output == T){

            message("Saving line lsmeans  and stas to folder in csv format")
            nm = paste(file.name, "_MEANS.csv", sep = "")
            nm.sts = paste(file.name, "_STATS.csv", sep = "")
            line.means = cbind(file.name, line.means)
            utils::write.csv(line.means, file = nm)
            utils::write.csv(sts, file = nm.sts)
        }

        means.sts.list = list(medias = as.data.frame(line.means), sts = sts)
        return(means.sts.list)
        print(sts)
    }



}



