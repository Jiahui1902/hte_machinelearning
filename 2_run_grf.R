
# begin by clearing all variables in current global environment 
# and prior sessions
rm(list=ls())

# part 1: administration commands -------------------------------------------------
# load necessary packages prior to run these codes below
pkgs <- c("haven","foreign","readstata13","data.table",
          "magrittr","stringr","devtools","rpart",
          "rpart.plot","rpart.utils","tidyverse","stargazer",
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

set.seed(4100)

# part 2: load data and run models -------------------------------------------------
# load data in R
edurose_mediation_20181126 <- read.dta13(
  "machinelearning/edurose_mediation_20181126.dta",convert.factors = F) %>%
  # transfer example dataset into data frame, which is the computation
  # unit used in R 
  as.data.frame

# The linear terms in this education dataset
# extract the covariantes to make things clearer when writing the
# formula in constructing the trees
linear_terms <- c("male","black","hisp","i_daded",
                  "i_momed","i_parinc","i_daduwhcol",
                  "i_intact","i_sibsz","i_rural",
                  "i_south","i_abil","i_hsprog",
                  "i_eduexp","i_eduasp","i_freduasp",
                  "i_rotter","i_delinq","i_schdisadv",
                  "i_mar18","i_parent18","good")




# If processed from Stata first, load the data file generated from Stata
# as well before running the codes
statamodel <- read.dta13("model_specification.dta")%>%
  as.data.frame

# define an empty list to save all the results
# this will make data export easier and reduce the number of variables 
# in the global environment
results <- list() # export results in rds file

################### model part ###################
ps_indicator <- 'propsc_com25_rf' # the propensity score from the 
outcomevariable <- 'lowwaprop'
treatment_indicator <- 'compcoll25'
# remove the observations where the covariates or propensity scores are missing
data.work <- edurose_mediation_20181126[,c(linear_terms,ps_indicator,
                                           'lowwaprop','compcoll25',
                                           "propcc3c","parinc3c",
                                           "momed3c","abil3c","race")]%>%
  na.omit

# decompose the matrix 
# covaraites matrix
X <- data.work[,c(linear_terms,ps_indicator)]
# outcome variable
Y <- data.work$lowwaprop
# treatment 
W <- data.work$compcoll25

# train causal forest model on the training set
tau.forest <- causal_forest(X, Y, W, num.trees = 4000)

# check the split point
lapply(tau.forest$`_split_vars`,function(xxx){
  unique(xxx)%>%length
})%>%unlist%>%table

##############################################
# Estimate treatment effects for the training data using out-of-bag prediction.
##############################################
tau.hat.oob = predict(tau.forest)$predictions

##############################################
# Estimate the conditional average treatment effect on the full sample (CATE).
estimates <- average_treatment_effect(tau.forest, target.sample = "all")

# Estimate the conditional average treatment effect on the treated sample (CATT).
average_treatment_effect(tau.forest, target.sample = "treated")
average_treatment_effect(tau.forest, target.sample = "all")

###########################
# best linear projection (average local treatment effect)
###########################

best_linear_projection(tau.forest,X)

best_linear_projection(tau.forest,X[,'propsc_com25_rf'])

# the results should be got:
# Best linear projection of the conditional average treatment effect.
# Confidence intervals are cluster- and heteroskedasticity-robust (HC3):
#   
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.216861   0.032754 -6.6209    4e-11 ***
# propsc_com25_rf 0.259032   0.091759  2.8230  0.00478 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

data.work$`estimated CATE` <- tau.hat.oob
vars <- c("propcc3c","parinc3c",
          "momed3c","abil3c","race")
############ extra table ##########
# identify the covariate stratified subgroups
rowsloc <- edurose_mediation_20181126%>%  # GRF does not allow missing values
  select_(.dots = as.list(c(outcomevariable,linear_terms,ps_indicator,treatment_indicator))) %>%
  na.omit() %>%
  rownames %>%
  as.numeric
table_extra <- do.call(rbind,lapply(1:length(vars),function(i){
  # get the subgroup difined by two variables
  grps <- edurose_mediation_20181126%>%  # GRF does not allow missing values
    .[rowsloc,] %>%
    mutate(group_id = group_indices_(., vars[i]))
  # get the estimation in each small square in the subgroup
  do.call(rbind,lapply(grps$group_id%>%unique%>%sort,function(xxx){
    subgroup <- grps$group_id==xxx
    x <- average_treatment_effect(tau.forest,
                                  target.sample = "all",subset = subgroup)
    return(c(x,Obs = sum(subgroup)))
  }))
})
)

# the results:
cbind(strata = paste0(rep(vars,each=3),
                      rep(1:3,5)),
      table_extra%>%as.data.frame)

