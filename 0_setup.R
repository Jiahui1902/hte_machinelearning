# This script is used for setting up the whole environment and
  # load necessary packages. In following scripts, it can be 
  # reused as the header in orther scripts.

# begin by clearing all variables in current global environment 
  # and prior sessions
rm(list=ls())
# Set 'working directory' to the folder where all the data and 
  # respective codes are located.

# to make output easier, record this working directory in one variable
wd <- getwd()

# part 1: administration commands -------------------------------------------------
# load necessary packages prior to run these codes below

pkgs <- c("haven","foreign","readstata13","data.table",
          "magrittr","stringr","devtools","rpart",
          "rpart.plot","latticeExtra","lattice","grid",
          "collapsibleTree","rpart.utils","dplyr",
          "Rlab","MASS","plot3D","mvtnorm","ggplot2",
          "Matrix",'Rlab','optmatch','reshape',
          'gtools','tidyverse',"htmlwidgets",
          'partykit','data.tree','grf','party','Matching')

# load packages: If the package is not sucessfully loaded, these codes
  # will intall this package first and then load it for another time.
installpackages <- lapply(pkgs,function(x){
if(x %in% rownames(installed.packages()) == FALSE) {install.packages(`x`)}
})
loadpackages <- lapply(pkgs,function(x){
library(`x`,character.only = T)
})
rm(list = c("loadpackages","installpackages"))

if (!require("causalTree")) devtools::install_github("susanathey/causalTree")
library("causalTree")


