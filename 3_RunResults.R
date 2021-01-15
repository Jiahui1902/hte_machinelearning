
# Part 0: set up and administration commands ----------------------------------------------------------
# begin by clearing all variables in current global environment and 
# prior sessions
rm(list=ls())

# use the header to set up working environment and load necessary 
# packages in R
# The working directory has been automatically set as the relavant folder
# to this .Rroj file.
source("machinelearning/1_setup.R")
# Export:1 - current working directory: wd
# 2 - directory for output data: outputdata
# 3 - directory for output graphs: gph

# check the output folders -- for the plots
# check if the folder for exporting results exist: if not, creat one
# (if so, this script will automatically skip this step)

if(!dir.exists("./outputdata")){
  dir.create("./outputdata")
}

if(!dir.exists(paste0("./outputdata/",Sys.Date(),"/"))){
  dir.create(paste0("./outputdata/",Sys.Date(),"/"))
}

if(!dir.exists("./gph")){
  dir.create("./gph")
}

if(!dir.exists(paste0("./gph/",Sys.Date(),"/"))){
  dir.create(paste0("./gph/",Sys.Date(),"/"))
}

# unlink(outputdata,recursive = T)

outputdata <- paste0("./outputdata/",Sys.Date(),"/")
gph <- paste0("./gph/",Sys.Date(),"/")

# part 1: administration commands -------------------------------------------------
# load data in R
edurose_mediation_20181126 <- read.dta13(
  # "rawdata/3c_edurose_mediation_20181126.dta") %>%
  "machinelearning/edurose_mediation_20181126.dta",convert.factors = F) %>%
  # transfer example dataset into data frame, which is the computation
  # unit used in R 
  as.data.frame


linear_terms <- c("male","black","hisp","i_daded",
                  "i_momed","i_parinc","i_daduwhcol",
                  "i_intact","i_sibsz","i_rural",
                  "i_south","i_abil","i_hsprog",
                  "i_eduexp","i_eduasp","i_freduasp",
                  "i_rotter","i_delinq","i_schdisadv",
                  "i_mar18","i_parent18","good")


# define an empty list to save all the results
# this will make data export easier and reduce the number of variables 
# in the global environment
results <- list() # export results in rds file

# Part II: run estimation algorithms -------------------------------------------------------------
# Construct causal tree with causal tree algorithms
# Load the functions used in constructing causal tree
source('machinelearning/1_CausalTree_Functions.R')

# define a dataframe including the potential outcome variables and the expected direction of larger treatment effect
femaleonlypool <-  list((1:nrow(edurose_mediation_20181126)),
                        (1:nrow(edurose_mediation_20181126))[which(edurose_mediation_20181126$male==0)]
)

{if(nrow(statamodel)==0){
  outcomevariables <- data.frame(vars=c("log_cum_income", # Cumulative wages
                                        "unempprop", # Proportion of time in unemployment
                                        "anywelf", # Receipt of welfare
                                        "vote06", # Voted in 2006
                                        "depress40", # Depressive symptoms
                                        "lowwaprop", # Prop. of time in low wage job 1990-2014
                                        "notmarriedprop",
                                        "eversinglenp",
                                        "transitions"
  ),minisize=20,crossvalidation=20,
  negative = as.logical(c(0,1,1,0,1,1,1,1,0)), # the expected direction of treatment effects
  femaleonly = as.numeric(c(0,0,0,0,0,0,1,1,0)+1)
  )
}
  else{
    outcomevariables <- statamodel[,c( unique( do.call(c,lapply(statamodel[1,],function(x){which(is.character(x))})) ),
                                       match("minsize",colnames(statamodel)),
                                       match("crossvalidation",colnames(statamodel)),
                                       match("negative",colnames(statamodel)),
                                       as.numeric( match("femaleonly",colnames(statamodel)) )
    )]
    outcomevariables$femaleonly <- outcomevariables$femaleonly+1
  }}



##################
# causal tree
##################
##################
# the original causal tree
##################
ps_indicator = 'propsc_com25_rf' # in order to get the covariates

set.seed(4382)
# original tree
hte_causalTree(outcomevariable = 'lowwaprop',minsize=20,crossvalidation = 40,negative = TRUE,
               data = edurose_mediation_20181126[femaleonlypool[[ as.numeric(outcomevariables[6,5]) ]], ],
               drawplot = TRUE,ps_indicator = 'propsc_com25_rf',no_indicater = '_minisize20',
               legend.x = 0.6,legend.y = 0.23)

# the one I need -------
results <- list() # export results in rds file

results_summary <- list() # export 
match_results <- list()


set.seed(4382) 
# adjusted by inversed propensity score weighting
hte_ipw(outcomevariable = 'lowwaprop',minsize=20,crossvalidation = 40,negative = TRUE,
               data = edurose_mediation_20181126[femaleonlypool[[ as.numeric(outcomevariables[6,5]) ]], ],
               drawplot = TRUE,ps_indicator = 'propsc_com25_rf',no_indicater = '_IPW',
               legend.x = 0.6,legend.y = 0.25)

# adjusted by grf
hte_forest(outcomevariable = "lowwaprop",
           minsize=40,crossvalidation = 40,
           negative = outcomevariables[7,4],
           data = edurose_mediation_20181126,
           drawplot = TRUE,
           ps_linear = "propsc_com25rflin",
           ps_indicator = 'propsc_com25_rf',
           treatment_indicator = "compcoll25",no_indicater = '_grf')


