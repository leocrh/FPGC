#' snptonumeric
#'@description  Recodes a single snp in hapmap format (double letters) to numeric format. For biallelic SNPs, homozygous are 0 and 2, and
#'heterozygotes are 1.
#'
#' @param snp A single SNP vector or dataframe in hapmap format, e.g., two letters for nucleotides without using characters for separation.
#'
#' @return Returns a SNP vector in numeric format
#'@export

snptonumeric <- function(snp = NULL) {

  f.snp = as.data.frame(prop.table(table(snp, useNA = "ifany")))
  f.snp  = f.snp[order(f.snp$Freq) ,]
  f.snp = as.data.frame(f.snp)

  f.snp
  f.snp = tidyr::separate(data = f.snp, col = snp, into = c("allele1", "allele2"), sep = 1,extra = "warn", remove = FALSE)
  f.snp$snp = as.character(f.snp$snp)

  f.snp$snp.num = NA
  f.snp = as.data.frame(f.snp)
  f.snp

  na.snps = f.snp[is.na(f.snp$snp) ,]

    het.snps = f.snp[!f.snp$allele1 == f.snp$allele2 & !is.na(f.snp$snp) ,]

  if(nrow(het.snps) == 1){
    het.snps$snp.num = 1
    het.snps
  }

  hm.snps = f.snp[f.snp$allele1 == f.snp$allele2 & !is.na(f.snp$snp),]
  hm.snps

  if(nrow(hm.snps) == 2){
    hm.snps = hm.snps[order(hm.snps$Freq) ,]
    hm.snps$snp.num = c(0, 2)
  }

  if(nrow(hm.snps) ==1) {
    hm.snps$snp.num = 2
  }

  f.snp = rbind(na.snps, hm.snps, het.snps)

  f.snp.num = f.snp[, names(f.snp) %in% c("snp", "snp.num")]

  snp = as.data.frame(snp)
  snp$id = rownames(snp)
  snp = snp[order(snp$id) ,]
  colnames(snp)[1] = "snp"
  snp.out = merge(x = snp, y = f.snp.num, by = "snp", sort = FALSE)
  snp.out = snp.out[order(snp.out$id) ,]
  rownames(snp.out) = rownames(snp)
  snp.out = snp.out[, 3]

  return(snp.out)
}

#' snpMatrixToNumeric
#' @description Recodes a biallelic SNP matrix in hapmap format to numeric format, e.g., 0, 1, 2. The function applies the snptonumeric function
#' over the columns of the SNP matrix.
#'
#' @param X A SNP matrix containing in its row names the id of the individuals and columns as SNPs.
#'
#' @return A numeric SNP matrix.
#' @seealso snptonumeric
#'@export
#'
snpMatrixToNumeric <- function(X = NULL) {

  X = as.matrix(X)
  X = X[order(rownames(X)) ,]

  X.num = apply(X, 2, snptonumeric)
  rownames(X.num) = rownames(X)
  return(X.num)

}
