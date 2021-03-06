% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gmatrix.R
\name{pcaStructure}
\alias{pcaStructure}
\title{PCA on G}
\usage{
pcaStructure(
  G,
  n.clusters = 3,
  plot.var = T,
  plot.pca = T,
  save.pca.plot = F,
  pca.plot.name = "PCA of G",
  ...
)
}
\arguments{
\item{G}{Is the relationship matrix}

\item{n.clusters}{The number of cluster to be solved over the kmeans algorithm}

\item{plot.var}{Logial to indicate if the scree plot is displayed}

\item{plot.pca}{Logical to indicate if the PCA plot is displayed with ggplot2}

\item{save.pca.plot}{Logical to indicate if the PCA plot will be saved}

\item{pca.plot.name}{Is the name of the file for the saved plot}

\item{...}{Arguments passed to the ggsave function}
}
\description{
Function to perform the PCA on the kinship matrix
}
