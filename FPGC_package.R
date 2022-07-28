# Step 0: Packages you will need
install.packages("devtools")

devtools::install_github("klutometis/roxygen")

library(roxygen2)
library("devtools")
library(usethis)
library(testthat)

# Step 1: Create your package directory
wd = ("C:/Users/LCRESPO/OneDrive - CIMMYT/Documents/FPGC")
setwd(wd)

create_package("FPGC")

#Step 2: Add functions

# Step 3: Add documentation

#Step 4: Process your documentation
setwd("./FPGC")
load_all()

use_package(package = "ggplot2",  type = "Imports", min_version = TRUE)
use_package(package = "lmerTest", type = "Imports", min_version = TRUE)
use_package(package = "lme4",  type = "Imports", min_version = TRUE)
use_package(package = "reshape2", type = "Imports", min_version = TRUE)

use_package(package = "nasapower", type = "Imports", min_version = TRUE)
use_package(package = "plyr", type = "Imports", min_version = TRUE)
use_package(package = "stringr", type = "Imports", min_version = TRUE)
use_package(package = "tidyr", type = "Imports", min_version = TRUE)
use_package(package = "doBy", type = "Imports", min_version = TRUE)
use_package(package = "zoo", type = "Imports", min_version = TRUE)

check()

?rmminorallele #example
?hapmaptonumeric_matrix
?snptonumeric
?snpMatrixToNumeric
?haptodoubl


document() #Create documentation

# Installing the package
setwd("..")
install("FPGC")

#Loading the package
devtools::install_github("leocrh/FPGC", force = T)
library(FPGC)
