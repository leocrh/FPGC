

#' Recodes single nucleotides to numeric format
#' @description This function takes one SNP marker with the single nucleotide hapmap format of snps to convert them to numeric format:
#' 0, 1, 2. Where 0 and 2 are the homozygous genotypes and 1 the herozygote genotypes. 0 is assigned to the major allele
#' and 2 to the minor allele.
#' @param gd.snp A single SNP vector in hapmap format
#'
#' @return Returns a SNP vector in numeric format
#' @export
#'
hapmaptonumeric <- function(gd.snp) {
  snp.table = as.data.frame(table(gd.snp)) #table of snp to calculate frequency
  snp.table$af = snp.table$Freq/sum(snp.table$Freq) #genotype frequency

  snp.table.hm = snp.table[snp.table$gd.snp %in% c("A", "G", "C", "T") ,] #keep homozygous genotypes
  snp.table.ht = snp.table[!snp.table$gd.snp %in% c("A", "G", "C", "T") ,] # Keep heterozygous gneotypes
  snp.table.hm = snp.table.hm[order(-snp.table.hm$af) ,]

  if(nrow(snp.table.hm) == 2){ # Piece of code when there are two homozygous genotypes

    if(nrow(snp.table.ht) == 1){ #Piece of code when there are heterozygotes and homozygous

      snp.table.ht$num = "1" # Setting heterozygotes as "1"
      snp.table.hm = snp.table.hm[order(-snp.table.hm$af) ,]
      snp.table.hm$num =  c("0", "2")
      snp.table = rbind(snp.table.hm, snp.table.ht)

      print(snp.table)

      t = setNames(snp.table$num, snp.table$gd.snp)
      t
      gd.num = gd.snp %>%
        stringr::str_replace_all(t)
      return(gd.num)

    }

    if(nrow(snp.table.ht) == 0) { #Piece of code when there are no heterozygotes
      snp.table.hm = snp.table.hm[order(-snp.table.hm$af) ,]
      snp.table.hm$num =  c("0", "2")
      t = setNames(snp.table.hm$num, snp.table.hm$gd.snp)
      t
      gd.num = gd.snp %>%
        stringr::str_replace_all(t)
      return(gd.num)
    }
  }

  if(nrow(snp.table.hm) == 1) { #Piece of code when there is one homozygous genotype
    snp.table.hm$num = "0"
    t = setNames(snp.table.hm$num, snp.table.hm$gd.snp)
    t
    gd.num = gd.snp %>%
      stringr::str_replace_all(t)
    return(gd.num)
  }

  if(nrow(snp.table.hm) == 0) { # Piece of code when there is only heterozygotes
    snp.table.hm$num = "1"
    t = setNames(snp.table.hm$num, snp.table.hm$gd.snp)
    t
    gd.num = gd.snp %>%
      stringr::str_replace_all(t)
    return(gd.num)
  }



}

####
#### Multple snp appliction of hapmaptonumeric
####

#' hapmaptonumeric_matrix
#' @description This function applyes the hapmaptonumeric function on a dataframe of SNP in hapmap format to recode as numeric: 0, 1, 2.
#' The marker data requires that the SNPs are in columns and the individual are in rows.
#'
#'
#' @param X SNP matrix in hapmapformat with single nucleotides
#' @param snp.start.col An index number where the SNP dat starts in X
#'
#' @return A dataframe of SNPs in numeric format
#' @export
#'
hapmaptonumeric_matrix = function(X = NULL, snp.start.col = 2) {
  X = X[, snp.start.col:ncol(X)]
  Xnum = apply(X = X, MARGIN = 2, hapmaptonumeric)
  rownames(Xnum) = rownames(X)
  Xnum = as.data.frame(Xnum)
  return(Xnum)
}
