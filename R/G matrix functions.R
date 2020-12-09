
##%######################################################%##
#                                                          #
####             Function to generate the G             ####
####       matrix from marker data. Three options       ####
####    are available: "cross.product" is simple the    ####
####         cross product of the marker data;          ####
####         "VanRaden" is the Van Raden method         ####
#### published in 2007; and "cross.product.average" is  ####
####          the cross product of the matrix           ####
####          average by the number of markers          ####
#                                                          #
##%######################################################%##

#' Generate G
#'
#' @description Function to generate the G matrix from marker data.
#' @param X Is a Marker matrix with marker Ids are in columns and Individuals are in rows
#' @param method One of  "cross.product", "cross.product.average" or "VanRaden". The default methods is "cross.product.average".
#' "cross.product is simply XX'.
#' "cross.product.average" is X centered and scaled and then XX' divided the number of markers
#' "VanRaden" is X centered and scaled by the allele frequency as in VanRaden 2007.
#' @export
#' @examples G = G.matrix(X = x, method = "VanRaden")


Gmatrix = function(X = x, method = "cross.product.average") {
  if(method == "cross.product") {
    G1=tcrossprod(as.matrix(X))
    G1
  }

  # method 2 is the Van Raden 2007 method
  else if(method == "VanRaden") {
    phat=colMeans(X, na.rm = T)/2
    MAF=ifelse(phat<0.5,phat,1-phat)
    MAF = as.data.frame(MAF)
    X2=scale(X,center=TRUE, scale=FALSE)
    k=2*sum(MAF*(1-MAF))
    G2=tcrossprod(X2)/k
    G2
  }

  else if(method == "cross.product.average") {
    X3=scale(X,center=TRUE,scale=TRUE)
    G3=tcrossprod(X3)/ncol(X3)
    G3

  }

}

#

##%######################################################%##
#                                                          #
####     Function to perform the PCA on the kinship     ####
####                  matrix. G is the                  ####
####      relationship matrix; n.clusters indicate      ####
####      the number of clusters to be solved over      ####
####           the kmeans algorithm; plot.var           ####
####              is logical to determine               ####
####      if the scree plot is output; plot.pca is      ####
####      logical to determine if the pca plot is       ####
####         output with ggplot                         ####
#                                                          #
##%######################################################%##

#' PCA on G
#'
#' @description Function to perform the PCA on the kinship matrix
#' @param G Is the relationship matrix
#' @param n.clusters The number of cluster to be solved over the kmeans algorithm
#' @param plot.var Logial to indicate if the scree plot is displayed
#' @param plot.pca Logical to indicate if the PCA plot is displayed with ggplot2
#' @param save.pca.plot Logical to indicate if the PCA plot will be saved
#' @param pca.plot.name Is the name of the file for the saved plot
#' @param ... Arguments passed to the ggsave function
#' @export

pcaStructure = function(G, n.clusters = 3, plot.var = T, plot.pca = T,  save.pca.plot = F,
                         pca.plot.name = "PCA of G", ...) {
  q.est = prcomp(G)
  message("SUMMARY OF THE PCA ON G")
  print(summary(q.est))
  pcs = summary(q.est)

  pc1 = pcs$importance[2,1] * 100
  pc1 = round(pc1, 2)
  pc2 = pcs$importance[2,2] * 100
  pc2 = round(pc2, 2)

  message("PROPORTION OF VARIANCE EXPLAINED BY THE FIRST TWO COMPONENTS")
  print(pc1)
  print(pc2)

  if(plot.var == T){
  plot(q.est, type = "l")}

  q.est.df = data.frame(q.est$x)
  fit1 <- kmeans(q.est.df, n.clusters, nstart = 100)
  q.est.df$cluster = as.factor(fit1$cluster)


  library(ggplot2)
  pca.plot = ggplot2::ggplot(q.est.df, aes(x = PC1, y = PC2, group = cluster, col = cluster)) +
  geom_hline(yintercept = 0, size = 0.1, linetype = "dashed") +
  geom_vline(xintercept = 0, size = 0.1, linetype = "dashed") +
  geom_point(size=3,alpha=0.8)+
  xlab(label = paste("PC1 Proportion of variance", pc1, sep = " ")) +
  ylab(label = paste("PC2 Proportion of variance", pc2, sep = " ")) +
  scale_y_continuous(breaks = seq(from = -20, to = 20, by = 2))+
  scale_x_continuous(breaks = seq(from = -20, to = 20, by = 2))+
  guides(colour = guide_legend("Group"))+
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "bottom")

  if(plot.pca == T) {print(pca.plot)}

  if(save.pca.plot == T) {
    ggsave(filename = paste(pca.plot.name, ".tiff", sep = ""), plot = pca.plot, path = getwd(),
          device = "tiff", ...)
  }

  return(q.est.df)

  }

