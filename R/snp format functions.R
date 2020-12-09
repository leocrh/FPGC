############################################################
#                                                          #
#       Function to convert HAPMAP format to double        #
# alleles with no separation snp.col indicates the column  #
#  number where the snps start is.gid.in.col indicates if  #
#      the genotypes are in columns to then transpose      #
#               the output   and generate a                #
#         matrix of genotypes in rows and snps in          #
#    columns snp.names is to indicate the column where     #
#        the snp names are if is.gid.in.col is TRUE        #
#                                                          #
############################################################

#' Function to convert HAPMAP format to double alleles with no separation, i.e, from A to AA, C to CC, R to AG, Y to CT, etc.,
#' @param df Is a data frame that contains SNP data in IUPAC format
#' @param snp.col Is the index where the SNP data starts in the data frame. Set to NULL if the SNPs are in rows
#' @param is.gid.in.col Is a Logical to indicate if the individuals are in columns (TRUE) or in rows (FALSE). Default value is FALSE
#' @param snp.names Is an index of the column where the SNP ids are if is.gid.in.col = TRUE

hap.to.doubl = function(df = df,
                        snp.col=NULL,
                        is.gid.in.col = FALSE,
                        snp.names = NULL) {
  df.sub = df[, snp.col:ncol(df)]
  df.sub = apply(df.sub, 2, as.character)
  df.sub[df.sub == "A"] = "AA"
  df.sub[df.sub == "G"] = "GG"
  df.sub[df.sub == "C"] = "CC"
  df.sub[df.sub == "T"] = "TT"

  df.sub[df.sub == "R"] = "AG"
  df.sub[df.sub == "Y"] = "CT"
  df.sub[df.sub == "S"] = "GC"
  df.sub[df.sub == "W"] = "AT"
  df.sub[df.sub == "K"] = "GT"
  df.sub[df.sub == "M"] = "AC"

  df.sub[df.sub == "N"] = NA

  df = cbind(df[,1:snp.col-1], df.sub)
  df[, snp.col:ncol(df)] = apply(df[, snp.col:ncol(df)], 2, as.factor)
  message("DIMENTIONS OF MARKER MATRIX")
  print(dim(df))
  df

  if(is.gid.in.col == TRUE) {
    df = t(df)
    colnames(df) = df[snp.names,]
    df = df[2:nrow(df),]
    message("DIMENTIONS OF MARKER MATRIX")
    print(dim(df))
    return(df)
  }


}



##%######################################################%##
#                                                          #
####      This function converts the hapmap format      ####
####         to numeric by susing the SNPassoc          ####
####          package. df is a dataframe where          ####
####         the markers are in columns and the         ####
####            gids are in rows. g.names is            ####
####            a vector of gids, taken from            ####
####          df. snp.col is the column number          ####
####         where the marker data starts in df         ####
#                                                          #
##%######################################################%##

#' This function converts the hapmap formatto numeric by susing the SNPassoc package.
#' @param df Is a dataframe where the markers are in columns and the gids are in rows
#' @param g.names Is a vector of genotype ids (gids) taken from df
#' @param snp.col Is the index where the marker data starts in the df


library(SNPassoc)

hpm.to.numeric.add = function(df = d, g.names = NULL, snp.col = NULL){
  df = SNPassoc::setupSNP(data = df, colSNPs = snp.col:ncol(df), sep = "")
  df.num = apply(df, 2, additive)
  row.names(df.num) = g.names
  message("Dimentions and Head of the numeric marker data")
  print(dim(df.num))
  print(df.num[1:5, 1:10])
  return(df.num)
}




############################################################
#                                                          #
#            This function removes the missing             #
#           marker data at certain threshold and           #
#                plots  a histogram of the                 #
#             missing data if plot.missing = T             #
#                                                          #
############################################################

#Threshold is a value in %

#' This function removes the missing marker data at certain threshold and plots  a histogram of the missing data
#' @param X Is a marker matrix with markers in columns and gids in rows
#' @param threshold Is a threshold for the % of the missing values in X
#' @param plot.missing Is a Logical to indicate if the the histogram of missing values is plotted. Default is False

rm.missing.snp = function(X = x, threshold = 20, plot.missing = F) {
 if(sum(is.na(X)) == 0){
   stop("There are no missing values in X, have they been imputed?")}

  per.missing = (colSums(is.na(X))/nrow(X))*100 #computes % of missing marker data
    missing.t = per.missing[per.missing<threshold] # keeps marker data with less than missing info
  X = X[, names(X) %in% names(missing.t)] # Subsets X to marker data < missing info
  print(X[1:5,1:5])
  X = as.data.frame(X)

  if(plot.missing == T) {
    hist(per.missing, col = "black", border = "white",
         main = "Histogram of missing markers",
         xlab = "Percentage of missing data", ylab = "Number of markers")
    X = as.data.frame(X)
    }
  print(dim(X))
  return(X)
}


#' This function is to remove the markers that have a determined proportion of heterozygotes
#' @param X Is a marker matrix with markers in columns and gids in rows
#' @param threshold Is a threshold for the % of the heterozygotes in X
#' @param plot.het Logical to indicate if the histogram of heterozygocity is displayed. Default is FALSE

rm.het.snp = function(X = x, threshold = 10, plot.het = F) {
  per.het = apply(gd.num, 2, function(x) {((sum(x[x==1]))/length(x))*100})
  het.keep = per.het[per.het < 10]
  X = as.data.frame(X)
  X = X[, names(X) %in% names(het.keep)]
  print(X[1:5,1:5])
  X = as.matrix(X)
  message("The new dimenssions of X are:")
  print(dim(X))

  if(plot.het == T) {
    hist(het.keep, col = "black", border = "white",
         main = "Histogram of heterozygocity",
         xlab = "Percentage of heterozygocity", ylab = "Number of markers")
    X = as.data.frame(X)
    }

    return(X)

}




############################################################
#                                                          #
#             This function removes markers at             #
#           certain minor allele freq treshold.            #
# If plot.minor.allele = T it plots a histogram of the MAF #
#                                                          #
############################################################

#' This function removes markers at certain minor allele freq treshold
#' @param X Is a marker matrix with markers in columns and gids in rows
#' @param minor.threshold Is a threshold for the minor allele frequency (MAF)
#' @param plot.minor.allele Logial to indicate if the histogram of the MAF is displayed

rm.minor.allele = function(X = x, minor.threshold = 0.05, plot.minor.allele = F) {
  phat=colMeans(X, na.rm = T)/2
  MAF=ifelse(phat<0.5, phat, 1-phat)
  MAF = as.data.frame(MAF)
  phat.maf = as.data.frame(subset(MAF, MAF > minor.threshold))
  colnames(phat.maf) = "phat"

  X.maf = X[,rownames(phat.maf)]
  print(X.maf[1:7,1:7])
  message("New dimensions of snp matrix")
  print(dim(X.maf))

  if(plot.minor.allele == T){
    hist(phat.maf$phat, main="Histogram of allele frequency",
         col = "black", border = "white",
         xlab = "Allele frequency", ylab = "Number of markers")
  }
  return(X.maf)
}


############################################################
#                                                          #
#       This two functions are to impute the missing       #
#       data in the marker matrix and keep a matrix        #
#     that is the result of the imputations. if there      #
#               are monomophic markers after               #
#             the imputation, they are removed             #
#                                                          #
############################################################


#Utility routines
#NOTE: This routine will give the list only with monomorphic markers with missing values,
#if you have monomorphic markers without missing values you have to remove it manually

#' Impute is a function to perfomr marker imputation and provde a list of monomorphic markers with missing values.
#' Monomorphic markers witout missing values have to be removed manually.
#' @param X Is a numeric marker matrix
#' @examples
#' X.impute = Impute(X)


Impute=function(X)
{
  monomorphic=numeric()
  for(i in 1:ncol(X))
  {
    cat('Imputing Marker ',i,'\n')
    {
      if(length(as.numeric(table(X[,i])))==1)
      {
        monomorphic=c(monomorphic,i)
      }else{
        tmp=table(X[,i])
        x=as.numeric(names(tmp))
        X[which(is.na(X[,i])),i]=sample(x=x,size=sum(is.na(X[,i])),replace=TRUE,prob=tmp/sum(tmp))
      }
    }
  }
  return(list(X=X,monomorphic=monomorphic))

}

#' this function takes the output of the Impute function to keep the adequate data matrix. Run after Impute
#' @param X the output of the Impute function
#' @param out is the out element in the list from the Impute function

XImp = function(X, out=out) {
  if (length(out$monomorphic) == 0){
    X = as.data.frame(out$X)
  }

  else {
    X = as.data.frame(out$X[,-out$monomorphic])
  }

}


