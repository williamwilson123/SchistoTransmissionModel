## create package containing C++ model

# load required packages 
library(Rcpp)
library(RcppNumerical)

## CREATE PACKAGE FOR TRANSMISSION MODEL ---------------------------------------
Rcpp.package.skeleton(name="SchistoTransmissionModel", cpp_files = "src/schistomod.cpp")

## Here, manually edit the DESCRIPTION file inside the package folder, as follows:

# Package: SchistoTransmissionModel
# Type: Package
# Title: A schistosomiasis transmission and immunity model
# Version: 1.1   ***incremement when update package***
# Date: 2024-02-05
# Author: Gregory Milne
# Maintainer: Gregory Milne <gmilne@rvc.ac.uk>
#   Description: A schistosomiasis transmission and immunity model.
# License: GPL (>= 2)
# Imports: Rcpp (>= 1.0.11)
# LinkingTo: Rcpp,RcppEigen,RcppNumerical
# Depends: Rcpp,RcppEigen,RcppNumerical
# Encoding: UTF-8
# VignetteBuilder: knitr
# VignetteIndexEntry: 'SchistoTransmissionModel'
# Suggests: knitr, rmarkdown

## NB: make sure to leave a blank line at the end

## Then add any auxillary R scripts into the package's 'R' folder
# e.g., add unpackingfunction.R which runs the model & unpacks nested into unnested output

## Edit NAMESPACE:
## Add:
# export(RunTransmissionModel)
## where RunTransmissionModel is the name of the R function

## Remove:
# exportPattern("^[[:alpha:]]+")

## Add vignette (.Rmd file) to folder called 'vignettes'

# Then add documentation & install package
setwd("SchistoTransmissionModel")  # Change working directory to your package directory

# ensure devtools up to date ( might be causing errors in document() )
# devtools::install_github("r-lib/devtools")

# check R tools PATH
# Sys.getenv("PATH")

#Add Rtools to PATH
# Sys.setenv(PATH = paste("C:/Rtools43/bin", Sys.getenv("PATH"), sep=";"))

devtools::document() #generate documentation
devtools::build(path=getwd())    #build the 'tarball' in current directory
devtools::install()  #install the package (if fails with '! System command 'Rcmd.exe' failed' error, restart R & try again)

# if vignettes not being found can install from tarball instead
install.packages("SchistoTransmissionModel_1.1.tar.gz", repos = NULL, type = "source")

# Load your package
library(SchistoTransmissionModel)

# Then you can use the Rcpp model via the R wrapper function like so: 
# SchistoTransmissionModel::RunTransmissionModel()

# Is the vignette available?
vignette(topic="SchistoTransmissionModel", package="SchistoTransmissionModel")

# if the main Rcpp function is still available after loading the package 
# (i.e., not just the R wrapper function), you can delete the NAMESPACE file and 
# generate one again using the below function:
# roxygen2::roxygenize()
