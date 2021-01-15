# Functions for running causal tree algorithms


# Function 1 --------------------------------------------------------------
# Ordinary causal tree: use the difference between streatment and 
# control group as the estimator
hte_causalTree <- function(outcomevariable, 
                           # the name of outcome variabls we are interested in
                           minsize=20, 
                           # minimum number of treated observations, 
                           # control observations in a leaf
                           crossvalidation = 20, 
                           # number of cross-validations to do
                           data = edurose_mediation_20181126, 
                           # can be changed, and the defaul one defined here 
                           # is edurose_mediation_20181126, the education dataset we are
                           # working on
                           treatment_indicator = 'compcoll25', # treatment variable
                           ps_indicator = 'propsc_com25', # propensity scores 
                           covariates = c(linear_terms,ps_indicator),
                           negative = FALSE, 
                           # can be changed, specify the expected direction 
                           # of the treatment effects
                           drawplot = TRUE, 
                           # export the graph of tree structure if true
                           no_indicater="",
                           legend.x = 0.08,
                           legend.y = 0.25,
                           ...){
  # delete all missing values which is required in causal tree model 
  # and use it as the train set in machine learning model
  trainset <-  data[!is.na(data[,outcomevariable]),]
  
  # set up the formula used for constructing causal tree
  # export the formula for causal tree model: Y~X
  formula <- as.formula(paste(outcomevariable," ~ ", 
                              paste(covariates,collapse = '+'), collapse= "+"))
  
  # set up propensity score
  if(nchar(ps_indicator)>0){
    # contruct tree
    tree <- causalTree(formula, 
                       # specify the model, outcome variable ~ covariates
                       data = trainset, # specify the dataset to be used
                       treatment = trainset[,treatment_indicator], 
                       # specify the treatment variable, must be 0/1 indicator
                       split.Rule = "CT", 
                       # specify split rule; for causal tree, use "CT" 
                       # NOTE: there are four different splitting rules, 
                       # they are different in the cross-validation criteria used
                       # to determine the tree structure
                       # 1 - TOT
                       # 2 - CT
                       # 3 - fit
                       # 4 - tstats
                       # 5 - totD
                       # 6 - ctD
                       # 7 - fitD
                       # 8 - tstatsD
                       cv.option = "CT", # specifify cross validation method 
                       # and there are four different methods -- tot, ct, fit, tstats
                       # for causal tree, use "CT"
                       split.Honest = T, cv.Honest = T, split.Bucket = F, 
                       
                       xval = crossvalidation, 
                       # number of cross-validations to do and the default number is 20
                       cp = 0, 
                       propensity = trainset[,ps_indicator],
                       # specify the propensity score; if is not specified, it will use sum(treatment) / nobs as the propensity score
                       minsize = minsize # minimum number of treated observations, control observations in a leaf
                       # the default minimum size is 20, according to Jennie and Yu Xie's paper (Estimating Heterogeneous Treatment Effects with Observational Data, 2012)
    )}else{
      tree <- causalTree(formula, 
                         data = trainset, 
                         treatment = trainset[,treatment_indicator], 
                         split.Rule = "CT", 
                         cv.option = "CT",
                         split.Honest = T, cv.Honest = F, split.Bucket = F, 
                         xval = crossvalidation, 
                         cp = 0, 
                         minsize = minsize 
      )
    }
  
  # prune this tree model to avoid the overfitting issues
  # get the complexity parameter (cp) to be trimmed--the least important splits
  opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
  # recursively snipping off the least important tree based on the complexity parameter (cp)
  opfit <- prune(tree, opcp)
  
  # return the predicted heterogeneous treatment effect
  hte_effect <- opfit$frame$yval[opfit$where]
  
  # if drawplots is TRUE, make plots and export the plots
  if(drawplot==TRUE){
    ttable <<- data.frame()
    makeplots(negative=negative, opfit.=opfit,gph=gph,trainset,
              covariates,outcomevariable,data,ttable,no_indicater,legend.x,legend.y)
  }else{
    print(c('Drawplot = ', drawplot))
  }
  
  # In observational study:
  # export the estimation results by adding three new indicators:
  # *_hi variable is 0/1 variable indicating whether the treatment effect has excess the threshhold
  # *_predictedTE is the predicted treatment effect (but causal tree approach predicts the effects of a treatment theoretically)
  # *_rank is the ranking from highest to lowest relative to other units
  
  # merge all the variables into one single data frame
  # first put them in different columns, and add the ID: R0000100
  
  # In simulation study: just export the estimation results for treatment effects
  if(identical(trainset,edurose_mediation_20181126)){
    output <- cbind(out_hi=0,
                    out_predict=hte_effect,
                    R0000100 = data$R0000100) %>%
      # change it as a data table and rerank the predicted heterogeneous treatment effect
      as.data.table %>%
      # rerank
      .[order(hte_effect,decreasing = TRUE)] %>%
      # add the rankings relative to other observations
      .[,paste0(outcomevariable,"_rank"):=1:nrow(.)]
    # rename the column names for this new data frame
    colnames(output)[1:2] <- paste0(outcomevariable,c("_hi","_predictedTE"))
    # make sure it is a data frame
    output <- as.data.frame(output)
    # define the variable: *_hi which is a 0/1 variable indicating whether the treatment effect has excess the threshhold
    # And the threshhold is set to be the top 20% highest treatment effect
    output[,paste0(outcomevariable,"_hi")] <- {if(negative==FALSE){
      ifelse(output[,paste0(outcomevariable,"_predictedTE")]>=quantile(unique(output[,paste0(outcomevariable,"_predictedTE")]),c(1,2,3,4,5)/5)[4],1,0)
    }else{
      ifelse(output[,paste0(outcomevariable,"_predictedTE")]<=quantile(unique(output[,paste0(outcomevariable,"_predictedTE")]),c(1,2,3,4,5)/5)[1],1,0)
    }}
  }else{
    output <- cbind(hte_effect)%>%
      `colnames<-`(paste0(outcomevariable,"_predictedTE"))%>%
      as.data.frame
  }
  
  
  # put the results in the result list
  # NOTE: "<<-" means this result is forced to be saved in the global environment
  # and this step is set to make sure this function could be used for many different outcome variables
  # we do not need to rewrite the function when a new outcome variable comes in
  results[[paste0(outcomevariable,"_causalTree_",length(results)+1)]] <<- output
  results[[paste0(outcomevariable,"_causalTree_Tree_",length(results)+1)]] <<- opfit
}


# Making plots function ---------------------------------------------------
# plot the tree model (heatmap)
makeplots <- function(negative,opfit.=opfit,gph,trainset,covariates,outcomevariable,data.=data,
                      hte_effect_setup,no_indicater="",legend.x=0.8,legend.y=0.25,...
){
  # opfit. <- opfit.
  # outcomevariable <- outcomevariable
  if(nrow(opfit.$frame)>1){
    pdf(paste0(gph,outcomevariable,"_causaltree",no_indicater,".pdf")) # the name could be changed
    # define the color palette to be used in this head map
    # rbPal <- colorRampPalette(c('blue','red')) # the color palette could be changed
    # rbPal <- colorRampPalette(c('gold','royal blue')) # the color palette could be changed
    rbPal <- colorRampPalette(c('golden rod','blue')) # the color palette could be changed
    # get the color for each of the leaves
    {if(negative==FALSE){
      od <- 1:length(opfit.$frame$yval)
    }else{
      od <- length(opfit.$frame$yval):1 
    }}
    
    cols <- rbPal(length(opfit.$frame$yval))[od]
    # change variables to its labels
    if(identical(colnames(data.),colnames(edurose_mediation_20181126))){
      opfit.$frame$var <- as.character(opfit.$frame$var)
      opfit.$frame$var[opfit.$frame$var%in%colnames(data.)] <-
        str_replace(attr(data.,"var.labels")[na.omit(match(opfit.$frame$var,colnames(edurose_mediation_20181126)))%>%as.numeric],"imputed ","")%>%
        str_replace(.,"Estimated ","")
    }else{
      # opfit.$frame$var <- as.character(opfit.$frame$var)
      opfit.$frame$var <- as.character(opfit.$frame$var)
      opfit.$frame$var[opfit.$frame$var%in%colnames(data.)] <-
        str_replace(attr(data.,"var.labels")[na.omit(match(opfit.$frame$var,colnames(edurose_mediation_20181126)))%>%as.numeric],"imputed ","")%>%
        str_replace(.,"Estimated ","")
    }
    
    # make plots
    if(length(ttable)>0){
      prp(opfit., # tree model estimated from causal tree
          type = 1, # could be delted. draw a split label at each split and a node label at each leaf
          # extra =100, # display extra information at the nodes, could be deleted
          nn = FALSE, # could be deleted. display the node numbers and the default value is FALSE
          box.palette = cols, # specify the colors in each node
          suffix = paste0(ttable$star,'\n(',ttable$se%>%round(.,3),
                          ')\n',ttable$SampleSize%>%round(.,1),"%"), # specify the label in each nodes
          yesno = 2, # specify if "yes" and "no" are shown
          yes.text = "y", # specify the text for "yes"
          no.text = "n", # specify the text for "no"
          trace = TRUE,
          varlen = 0,
          col = "white",
          # branch.col = cols, # specify branch colors
          # split.col = cols, # specify the color of the split label text
          # nn.col = cols,# specify the color of the node numbers
          main=paste0("Heterogeneous Treatment Effects: ",
                      attr(edurose_mediation_20181126,"var.labels")[na.omit(match(outcomevariable,colnames(edurose_mediation_20181126)))%>%as.numeric]) # specify the title displayed at the top
      )
    }else{
      prp(opfit., # tree model estimated from causal tree
          type = 1, # could be delted. draw a split label at each split and a node label at each leaf
          extra =100, # display extra information at the nodes, could be deleted
          nn = FALSE, # could be deleted. display the node numbers and the default value is FALSE
          box.palette = cols, # specify the colors in each node
          yesno = 2, # specify if "yes" and "no" are shown
          yes.text = "y", # specify the text for "yes"
          no.text = "n", # specify the text for "no"
          trace = TRUE,
          varlen = 0,
          col = "white",
          # branch.col = cols, # specify branch colors
          # split.col = cols, # specify the color of the split label text
          # nn.col = cols,# specify the color of the node numbers
          main=paste0("Heterogeneous Treatment Effects: ",
                      attr(edurose_mediation_20181126,"var.labels")[na.omit(match(outcomevariable,colnames(edurose_mediation_20181126)))%>%as.numeric]) # specify the title displayed at the top
      )
    }
    # legend and notes
    mtext(c("\nText in Squares: HTE & sample size (%); Color of Squares: Blue:largest treatment effects & Yellow:smallest treatment effects; Number in Parentheses: Standard Error"),
          cex=0.5)
    legend(legend.x,legend.y, legend = c("HTE","sample size (%)",
                                 "largest treatment effects",
                                 "smallest treatment effects"),
           title = "Numbers/Colors:",cex = .6,box.lty=0,#pch = c(0,0,0,0),
           fill=c("white","white","blue","golden rod"))
    
    # alternatively, use the following command to make the plot
    # rpart.plot(opfit.,roundint = FALSE,type = 4,col = cols,main = paste0("Heterogeneous Treatment Effects: ",outcomevariable))
    dev.off()
    
    results[[paste0(outcomevariable,"_color_",length(results)+1)]] <<- opfit.
    cols <- cols[order(order(opfit.$frame$yval))]
    results[[paste0(outcomevariable,"_cols_",length(results)+1)]] <<- cols
    
    
    
    #   # animated results
    # pdf(paste0(gph,outcomevariable,"_animated_tree.pdf"),width = 13.33,height = 7.5)
    #   # cols <- cols[order(order(opfit.$frame$yval))]
    #   for(iframe in 1:nrow(opfit.$frame)) {
    #   # for(iframe in sort(rownames(opfit.$frame)%>%as.numeric)) {
    #     
    #     cols_new <- ifelse(1:nrow(opfit.$frame) <= iframe, cols, "white")
    #     cols_text <- ifelse(1:nrow(opfit.$frame) <= iframe, "black", "white")
    #     cols_text_yesno <- ifelse(1:nrow(opfit.$frame) < iframe, "black", "white")
    #     cols_text_yesno <- ifelse(opfit.$frame$var=="<leaf>", "white", cols_text_yesno)
    #     
    #     # cols_new <- ifelse(rownames(opfit.$frame)%>%as.numeric <= iframe, cols, "white")
    #     # cols_text <- ifelse(rownames(opfit.$frame)%>%as.numeric <= iframe, "black", "white")
    #     # cols_text_yesno <- ifelse(rownames(opfit.$frame)%>%as.numeric < iframe, "black", "white")
    #     # cols_text_yesno <- ifelse(opfit.$frame$var=="<leaf>", "white", cols_text_yesno)
    #     
    #     # make plots
    #     if(length(ttable)>0){
    #     prp(opfit., # tree model estimated from causal tree
    #         type = 1, # could be delted. draw a split label at each split and a node label at each leaf
    #         # extra =100, # display extra information at the nodes, could be deleted
    #         # nn = FALSE, # could be deleted. display the node numbers and the default value is FALSE
    #         box.col = cols_new, # specify the colors in each node
    #         suffix = paste0(ttable$star,'\n(',ttable$se%>%round(.,3),
    #                         ')\n',ttable$SampleSize%>%round(.,1),"%","\n\nY\t\t\t\t\tN"), # specify the label in each nodes
    #         under.cex = 1.5,
    #         under.col = cols_text_yesno,
    #         yesno = 0, # specify if "yes" and "no" are shown
    #         yes.text = "y", # specify the text for "yes"
    #         no.text = "n", # specify the text for "no"
    #         trace = TRUE,
    #         varlen = 0,
    #         col = "white",
    #         branch.col = cols_text, # specify branch colors
    #         split.col = cols_text, # specify the color of the split label text
    #         nn.col = cols_text,# specify the color of the node numbers
    #         # split.border.col= cols_text,
    #         main=paste0("Heterogeneous Treatment Effects: ",
    #                     attr(edurose_mediation_20181126,"var.labels")[na.omit(match(outcomevariable,colnames(edurose_mediation_20181126)))%>%as.numeric])
    #         )
    #     }else{
    #       prp(opfit., # tree model estimated from causal tree
    #           type = 1, # could be delted. draw a split label at each split and a node label at each leaf
    #           extra =100, # display extra information at the nodes, could be deleted
    #           # nn = FALSE, # could be deleted. display the node numbers and the default value is FALSE
    #           box.col = cols_new, # specify the colors in each node
    #           under.cex = 1.5,
    #           under.col = cols_text_yesno,
    #           yesno = 0, # specify if "yes" and "no" are shown
    #           yes.text = "y", # specify the text for "yes"
    #           no.text = "n", # specify the text for "no"
    #           trace = TRUE,
    #           varlen = 0,
    #           col = "white",
    #           branch.col = cols_text, # specify branch colors
    #           split.col = cols_text, # specify the color of the split label text
    #           nn.col = cols_text,# specify the color of the node numbers
    #           # split.border.col= cols_text,
    #           main=paste0("Heterogeneous Treatment Effects: ",
    #                       attr(edurose_mediation_20181126,"var.labels")[na.omit(match(outcomevariable,colnames(edurose_mediation_20181126)))%>%as.numeric])
    #       )
    #     }
    #     # legend and notes
    #     mtext(c("\nText in Squares: HTE & sample size (%); Color of Squares: Red:largest treatment effects & Blue:smallest treatment effects; Number in Parentheses: Standard Error"),
    #           cex=0.5)
    #     legend(legend.x,legend.y, legend = c("HTE","sample size (%)",
    #                                  "largest treatment effects",
    #                                  "smallest treatment effects"),
    #            title = "Numbers/Colors:",cex = .6,box.lty=0,#pch = c(0,0,0,0),
    #            fill=c("white","white","red","blue"))
    #   }
    # dev.off()
    
    
    
  }else{
    print(outcomevariable)
  }
}


# Function 2 --------------------------------------------------------------
# set up environment for matching
options("optmatch_max_problem_size" = Inf)
# matching on leaves
hte_matchinginleaves <- function(outcomevariable, 
                                 # the name of outcome variabls we are interested in
                                 minsize=20, 
                                 # minimum number of treated observations, 
                                 # control observations in a leaf
                                 crossvalidation = 20, 
                                 # number of cross-validations to do
                                 data = edurose_mediation_20181126, 
                                 # can be changed, and the defaul one defined here 
                                 # is edurose_mediation_20181126, the education dataset we are
                                 # working on
                                 treatment_indicator = 'compcoll25', # treatment variable
                                 ps_indicator = 'propsc_com25', # propensity scores 
                                 ps_linear = 'propsc_com25lin',
                                 covariates = c(linear_terms,ps_indicator),
                                 negative = FALSE, 
                                 # can be changed, specify the expected direction 
                                 # of the treatment effects
                                 drawplot = TRUE, 
                                 # export the graph of tree structure if true
                                 con.num=1,
                                 # the number of control variables used in matching
                                 no_indicater="",
                                 legend.x = 0.08,
                                 legend.y = 0.25,
                                 ...){
  # delete all missing values which is required in causal tree model 
  # and use it as the train set in machine learning model
  trainset <-  data[!is.na(data[,outcomevariable]),]
  
  # set up the formula used for constructing causal tree
  # export the formula for causal tree model: Y~X
  if(nchar(ps_indicator)>0){
    covariates_ <- c(covariates,ps_indicator) # non-linear ps score
    formula <- as.formula(paste(outcomevariable," ~ ", 
                                paste(covariates,collapse = '+'), collapse= "+"))
    covariates <- c(covariates,ps_linear) # linear ps score
  }else{
    formula <- as.formula(paste(outcomevariable," ~ ", 
                                paste(covariates,collapse = '+'), collapse= "+"))
    covariates <- covariates
  }
  # set up propensity score
  if(nchar(ps_indicator)>0){
    # contruct tree
    tree <- causalTree(formula, 
                       # specify the model, outcome variable ~ covariates
                       data = trainset, # specify the dataset to be used
                       treatment = trainset[,treatment_indicator], 
                       # specify the treatment variable, must be 0/1 indicator
                       split.Rule = "CT", 
                       # specify split rule; for causal tree, use "CT" 
                       # NOTE: there are four different splitting rules, 
                       # they are different in the cross-validation criteria used
                       # to determine the tree structure
                       # 1 - TOT
                       # 2 - CT
                       # 3 - fit
                       # 4 - tstats
                       # 5 - totD
                       # 6 - ctD
                       # 7 - fitD
                       # 8 - tstatsD
                       cv.option = "CT", # specifify cross validation method 
                       # and there are four different methods -- tot, ct, fit, tstats
                       # for causal tree, use "CT"
                       split.Honest = T, cv.Honest = T, split.Bucket = F, 
                       
                       xval = crossvalidation, 
                       # number of cross-validations to do and the default number is 20
                       cp = 0, 
                       propensity = trainset[,ps_indicator],
                       # specify the propensity score; if is not specified, it will use sum(treatment) / nobs as the propensity score
                       minsize = minsize # minimum number of treated observations, control observations in a leaf
                       # the default minimum size is 20, according to Jennie and Yu Xie's paper (Estimating Heterogeneous Treatment Effects with Observational Data, 2012)
    )}else{
      tree <- causalTree(formula, 
                         data = trainset, 
                         treatment = trainset[,treatment_indicator], 
                         split.Rule = "CT", 
                         cv.option = "CT",
                         split.Honest = T, cv.Honest = F, split.Bucket = F, 
                         xval = crossvalidation, 
                         cp = 0, 
                         minsize = minsize 
      )
    }
  
  # prune this tree model to avoid the overfitting issues
  # get the complexity parameter (cp) to be trimmed--the least important splits
  opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
  # recursively snipping off the least important tree based on the complexity parameter (cp)
  opfit <- prune(tree, opcp)
  
  # matchin in leaves and return the predicted heterogeneous treatment effect
  prty <- as.party(opfit)
  opfit_tree <- as.Node(prty)
  
  hte_effect_setup <- list()
  # Matching Algorithms
  x_num <- 0
  opfit_tree$Do(function(node){
    # extract data
    match_data <- data[rownames(node$data)%>%as.numeric,]
    match_number <- con.num
    x_num <<- x_num+1
    hte_effect_help. <- matchinleaves(trainset=match_data, # endogeneous
                                      covariates=covariates,
                                      outcomevariable=outcomevariable,
                                      hte_effect_setup = hte_effect_setup,
                                      treatment_indicator = treatment_indicator,
                                      con.num = match_number, # the numbers of controls in doing matching
                                      stars_setup)
    print(hte_effect_help.)
    hte_effect_setup[[x_num]] <<- c(hte_effect_help.,
                                    nrow(match_data)/nrow(trainset)*100%>%round(.,1))
    # keep all related numbers in the environment of node
    node$predicted <- hte_effect_help.[1]
    node$pvalue <- hte_effect_help.[2]
    node$standarderror <- hte_effect_help.[3]
    node$samplesize <- hte_effect_help.[4]
  })
  opfit_tree<<-opfit_tree
  hte_effect <- opfit_tree$Get("predicted")%>%as.numeric
  opfit$frame$yval <- hte_effect
  
  # create a new variable indicating the estimated treatment effect for each unit
  hte_effect <- opfit$frame$yval[opfit$where]
  
  # statistics
  ttable <<- unlist(hte_effect_setup)%>%matrix(.,ncol = 4,byrow = TRUE)%>%
    as.data.frame%>%
    `colnames<-`(c("Estimator","pvalue","se","SampleSize"))
  st <- rep("",length(ttable$pvalue))
  st[ttable$pvalue<0.05] <- "*"
  st[ttable$pvalue<0.01] <- "**"
  st[ttable$pvalue<0.001] <- "***"
  ttable$star <<- st
  # ttable$star <<- ifelse(ttable$pvalue<0.1,"*","")
  ttable$SampleSize <- round(ttable$SampleSize,1)
  
  # If makeing plots, the values from the original tree should be 
  # adjusted to the value generated from matching methods
  # adj_effect <- table(hte_effect)%>%as.data.table
  # opfit$frame$yval[match(adj_effect$N,opfit$frame$n)] <- as.numeric(adj_effect$hte_effect)
  
  
  # if drawplots is TRUE, make plots and export the plots
  if(drawplot==TRUE){
    # makeplots(opfit,gph,trainset,covariates,outcomevariable)
    makeplots(negative=negative, opfit.=opfit,gph=gph,trainset,
              covariates,outcomevariable,data.=data,ttable,no_indicater,legend.x,legend.y)
  }else{
    print(c('Drawplot = ', drawplot))
  }
  
  # In observational study:
  # export the estimation results by adding three new indicators:
  # *_hi variable is 0/1 variable indicating whether the treatment effect has excess the threshhold
  # *_predictedTE is the predicted treatment effect (but causal tree approach predicts the effects of a treatment theoretically)
  # *_rank is the ranking from highest to lowest relative to other units
  
  # merge all the variables into one single data frame
  # first put them in different columns, and add the ID: R0000100
  
  # In simulation study: just export the estimation results for treatment effects
  if(identical(trainset,edurose_mediation_20181126)){
    output <- cbind(out_hi=0,
                    out_predict=hte_effect,
                    R0000100 = data$R0000100) %>%
      # change it as a data table and rerank the predicted heterogeneous treatment effect
      as.data.table %>%
      # rerank
      .[order(hte_effect,decreasing = TRUE)] %>%
      # add the rankings relative to other observations
      .[,paste0(outcomevariable,"_rank"):=1:nrow(.)]
    # rename the column names for this new data frame
    colnames(output)[1:2] <- paste0(outcomevariable,c("_hi","_predictedTE"))
    # make sure it is a data frame
    output <- as.data.frame(output)
    # define the variable: *_hi which is a 0/1 variable indicating whether the treatment effect has excess the threshhold
    # And the threshhold is set to be the top 20% highest treatment effect
    output[,paste0(outcomevariable,"_hi")] <- {if(negative==FALSE){
      ifelse(output[,paste0(outcomevariable,"_predictedTE")]>=
               quantile(unique(output[,paste0(outcomevariable,"_predictedTE")]),c(1,2,3,4,5)/5)[4],1,0)
    }else{
      ifelse(output[,paste0(outcomevariable,"_predictedTE")]<=
               quantile(unique(output[,paste0(outcomevariable,"_predictedTE")]),c(1,2,3,4,5)/5)[1],1,0)
    }}
  }else{
    output <- cbind(hte_effect)%>%
      `colnames<-`(paste0(outcomevariable,"_predictedTE"))%>%
      as.data.frame
  }
  
  # put the results in the result list
  # NOTE: "<<-" means this result is forced to be saved in the global environment
  # and this step is set to make sure this function could be used for many different outcome variables
  # we do not need to rewrite the function when a new outcome variable comes in
  results[[paste0(outcomevariable,"_causalTree_matching_",length(results)+1)]] <<- output
  results[[paste0(outcomevariable,"_causalTree_matching_Tree_",length(results)+1)]] <<- opfit
  results[[paste0(outcomevariable,"_causalTree_matching_ttable_",length(results)+1)]] <<- ttable
}


###########################################
# matching function: 1:1 matching, full matching, 1:4 matching
###########################################
matchinleaves <- function(trainset=match_data,
                          covariates=covariates,
                          outcomevariable=outcomevariable,
                          hte_effect_setup = hte_effect_setup,
                          treatment_indicator,
                          con.num=1, # the numbers of controls in doing matching
                          ...){
  # variables used in the function
  outcomevariable <- outcomevariable
  con.num <- con.num
  
  # make sure the matrix for predicting propensity score is of full rank
  covariates.match <- covariates[apply(trainset[,covariates],2,
                                       function(i) unique(i)%>%length>1)]
  # set up the matching rule
  # match.formula <- as.formula(paste(treatment_indicator,'~',
  # paste(covariates.match,collapse = '+')))
  match.formula <- as.formula(paste(treatment_indicator,'~',covariates[length(covariates)]))
  
  print(match.formula)
  ######
  # matching in the leaves
  ######
  # optimal matching in leaves
  if(con.num==1){
    print("con.num=1")
    match <- optmatch::pairmatch(match_on(match.formula,
                                          # choose methods
                                          # step 1: assign the varible that should be matched exactly, and here is the propensity score strata
                                          # matchit can not deal with data with missing values. So I used fill.NAs functions to non-
                                          # informatively fill in missing values in original data frame. And these data are prepared to
                                          # to do optimal match.
                                          # step 2: create treatment to control distance with match_on function
                                          # step 3: use fullmatch to complete this optimal matching
                                          # notes: exactMatch can be neglected if units are allowed to be matched across levels
                                          method = "optimal",data = fill.NAs(trainset),
                                          tol=0.001)
                                 # controls = n (Optional, for example, set n = 1) 
                                 # The number of controls used in matching process. But this constraints may not 
                                 # be feasible. If so, R will automatically change the constraints and an warning message 
                                 # will be issued. 
                                 ,data = fill.NAs(trainset))
  }else if(con.num==4){
    # delete rows that contains missing values
    matchdata_nomiss <-  trainset%>%  # MatchIt does not allow missing values
      select_(.dots = as.list(c(outcomevariable,covariates[length(covariates)],treatment_indicator))) %>%
      na.omit()
    
    # the default setting is matching without replacement 
    # (by setting replacement to TRUE to do matching with replacement)
    X  <- matchdata_nomiss[,covariates[length(covariates)]]#%>%cbind
    Y  <- matchdata_nomiss[,outcomevariable]
    Tr  <- matchdata_nomiss[,treatment_indicator]
    
    match  <- Matching::Match(Y=Y, Tr=Tr, X=X, estimand = "ATE",
                    M=con.num,Weight = 2,replace=T)
    
  }else{
    match <- optmatch::fullmatch(match.formula,
                                 data = fill.NAs(trainset))
  }
  
  match_results[[length(match_results)+1]] <<- match
  
  # get the results
  if(con.num!=4){
    # get the new dataset after matching
    data.match <- cbind(fill.NAs(trainset),match)
    # get the treatment effects
    
    average <- c(paste0('mean(',outcomevariable,')'))
    average_name <- c(outcomevariable)
    # print(head(data.match))
    hte_effect_help <- data.match[which(!is.na(data.match$match)),
                                  c(treatment_indicator,outcomevariable,'match')]%>%
      group_by_(.dots = c('match',treatment_indicator))%>%
      summarise_(.dots = setNames(average, average_name))%>%
      as.data.frame
    hte_effect_help <- reshape(hte_effect_help,idvar='match',timevar = treatment_indicator,direction = 'wide')
    
    # extract the control and treatment group
    treatment <- hte_effect_help[,paste0(outcomevariable,'.1')]
    control <- hte_effect_help[,paste0(outcomevariable,'.0')]
    
    hte_effect_help <- mean(hte_effect_help[,paste0(outcomevariable,'.1')]
                            -hte_effect_help[,paste0(outcomevariable,'.0')])
    
    star <- t.test(treatment-control,conf.level = .9)
    # pvalue
    pvalue <- star$p.value
    # standard error
    standarderror <<- sqrt( var(treatment-control,na.rm = T)/length(na.omit(treatment-control)) )
  }else {
    pvalue <- 2*(1-pnorm(abs(match$est/match$se)))
    standarderror <<- match$se
    hte_effect_help <- match$est
  }
  print(hte_effect_help)
  print(c('The number of units forget to consider is',length(hte_effect_setup%>%.[is.na(.)])))
  c(hte_effect_help,pvalue,standarderror)%>%round(.,4)
}



#####################################
# Inverse Propensity Score Weighting#
#####################################

hte_ipw <- function(outcomevariable, 
                                 # the name of outcome variabls we are interested in
                                 minsize=20, 
                                 # minimum number of treated observations, 
                                 # control observations in a leaf
                                 crossvalidation = 20, 
                                 # number of cross-validations to do
                                 data = edurose_mediation_20181126, 
                                 # can be changed, and the defaul one defined here 
                                 # is edurose_mediation_20181126, the education dataset we are
                                 # working on
                                 treatment_indicator = 'compcoll25', # treatment variable
                                 ps_indicator = 'propsc_com25', # propensity scores 
                                 ps_linear = 'propsc_com25lin',
                                 covariates = c(linear_terms,ps_indicator),
                                 negative = FALSE, 
                                 # can be changed, specify the expected direction 
                                 # of the treatment effects
                                 drawplot = TRUE, 
                                 # export the graph of tree structure if true
                                 con.num=1,
                                 # the number of control variables used in matching
                                 no_indicater="",
                                 legend.x = 0.08,
                                 legend.y = 0.25,
                                 ...){
  # delete all missing values which is required in causal tree model 
  # and use it as the train set in machine learning model
  trainset <-  data[!is.na(data[,outcomevariable]),]
  
  # set up the formula used for constructing causal tree
  # export the formula for causal tree model: Y~X
  if(nchar(ps_indicator)>0){
    covariates_ <- c(covariates,ps_indicator) # non-linear ps score
    formula <- as.formula(paste(outcomevariable," ~ ", 
                                paste(covariates,collapse = '+'), collapse= "+"))
    covariates <- c(covariates,ps_linear) # linear ps score
  }else{
    formula <- as.formula(paste(outcomevariable," ~ ", 
                                paste(covariates,collapse = '+'), collapse= "+"))
    covariates <- covariates
  }
  # set up propensity score
  if(nchar(ps_indicator)>0){
    # contruct tree
    tree <- causalTree(formula, 
                       # specify the model, outcome variable ~ covariates
                       data = trainset, # specify the dataset to be used
                       treatment = trainset[,treatment_indicator], 
                       # specify the treatment variable, must be 0/1 indicator
                       split.Rule = "CT", 
                       # specify split rule; for causal tree, use "CT" 
                       # NOTE: there are four different splitting rules, 
                       # they are different in the cross-validation criteria used
                       # to determine the tree structure
                       # 1 - TOT
                       # 2 - CT
                       # 3 - fit
                       # 4 - tstats
                       # 5 - totD
                       # 6 - ctD
                       # 7 - fitD
                       # 8 - tstatsD
                       cv.option = "CT", # specifify cross validation method 
                       # and there are four different methods -- tot, ct, fit, tstats
                       # for causal tree, use "CT"
                       split.Honest = T, cv.Honest = T, split.Bucket = F, 
                       
                       xval = crossvalidation, 
                       # number of cross-validations to do and the default number is 20
                       cp = 0, 
                       propensity = trainset[,ps_indicator],
                       # specify the propensity score; if is not specified, it will use sum(treatment) / nobs as the propensity score
                       minsize = minsize # minimum number of treated observations, control observations in a leaf
                       # the default minimum size is 20, according to Jennie and Yu Xie's paper (Estimating Heterogeneous Treatment Effects with Observational Data, 2012)
    )}else{
      tree <- causalTree(formula, 
                         data = trainset, 
                         treatment = trainset[,treatment_indicator], 
                         split.Rule = "CT", 
                         cv.option = "CT",
                         split.Honest = T, cv.Honest = F, split.Bucket = F, 
                         xval = crossvalidation, 
                         cp = 0, 
                         minsize = minsize 
      )
    }
  
  # prune this tree model to avoid the overfitting issues
  # get the complexity parameter (cp) to be trimmed--the least important splits
  opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
  # recursively snipping off the least important tree based on the complexity parameter (cp)
  opfit <- prune(tree, opcp)
  
  # matchin in leaves and return the predicted heterogeneous treatment effect
  prty <- as.party(opfit)
  opfit_tree <- as.Node(prty)
  
  hte_effect_setup <- list()
  # Matching Algorithms
  x_num <- 0
  opfit_tree$Do(function(node){
    x_num <<- x_num+1
    # extract data
    ipw_data <- data[(rownames(node$data)%>%as.numeric),]
    wt <- (ipw_data[,treatment_indicator] - ipw_data[,ps_indicator])/(ipw_data[,ps_indicator] * (1-ipw_data[,ps_indicator]))
    # fit the model
    fit <- lm(as.formula(paste0(outcomevariable,'~',treatment_indicator)),
              data = ipw_data,
              weights = abs(wt))
    print(summary(fit)$coef[2,])
    hte_effect_setup[[x_num]] <<- cbind(summary(fit)$coef[2,c(1,4,2)]%>%t,
                                    nrow(node$data)/nrow(trainset)*100%>%round(.,1))
    # keep all related numbers in the environment of node
    node$predicted <- summary(fit)$coef[2,1]
    node$pvalue <- summary(fit)$coef[2,4]
    node$standarderror <- summary(fit)$coef[2,2]
    node$samplesize <- nrow(node$data)
  })
  opfit_tree<<-opfit_tree
  hte_effect <- opfit_tree$Get("predicted")%>%as.numeric
  opfit$frame$yval <- hte_effect
  
  # create a new variable indicating the estimated treatment effect for each unit
  hte_effect <- opfit$frame$yval[opfit$where]
  
  # statistics
  ttable <<- unlist(hte_effect_setup)%>%matrix(.,ncol = 4,byrow = TRUE)%>%
    as.data.frame%>%
    `colnames<-`(c("Estimator","pvalue","se","SampleSize"))
  st <- rep("",length(ttable$pvalue))
  st[ttable$pvalue<0.05] <- "*"
  st[ttable$pvalue<0.01] <- "**"
  st[ttable$pvalue<0.001] <- "***"
  ttable$star <<- st
  # ttable$star <<- ifelse(ttable$pvalue<0.1,"*","")
  ttable$SampleSize <- round(ttable$SampleSize,1)
  
  # If makeing plots, the values from the original tree should be 
  # adjusted to the value generated from matching methods
  # adj_effect <- table(hte_effect)%>%as.data.table
  # opfit$frame$yval[match(adj_effect$N,opfit$frame$n)] <- as.numeric(adj_effect$hte_effect)
  
  
  # if drawplots is TRUE, make plots and export the plots
  if(drawplot==TRUE){
    # makeplots(opfit,gph,trainset,covariates,outcomevariable)
    makeplots(negative=negative, opfit.=opfit,gph=gph,trainset,
              covariates,outcomevariable,data.=data,ttable,no_indicater,legend.x,legend.y)
  }else{
    print(c('Drawplot = ', drawplot))
  }
  
  # In observational study:
  # export the estimation results by adding three new indicators:
  # *_hi variable is 0/1 variable indicating whether the treatment effect has excess the threshhold
  # *_predictedTE is the predicted treatment effect (but causal tree approach predicts the effects of a treatment theoretically)
  # *_rank is the ranking from highest to lowest relative to other units
  
  # merge all the variables into one single data frame
  # first put them in different columns, and add the ID: R0000100
  
  # In simulation study: just export the estimation results for treatment effects
  if(identical(trainset,edurose_mediation_20181126)){
    output <- cbind(out_hi=0,
                    out_predict=hte_effect,
                    R0000100 = data$R0000100) %>%
      # change it as a data table and rerank the predicted heterogeneous treatment effect
      as.data.table %>%
      # rerank
      .[order(hte_effect,decreasing = TRUE)] %>%
      # add the rankings relative to other observations
      .[,paste0(outcomevariable,"_rank"):=1:nrow(.)]
    # rename the column names for this new data frame
    colnames(output)[1:2] <- paste0(outcomevariable,c("_hi","_predictedTE"))
    # make sure it is a data frame
    output <- as.data.frame(output)
    # define the variable: *_hi which is a 0/1 variable indicating whether the treatment effect has excess the threshhold
    # And the threshhold is set to be the top 20% highest treatment effect
    output[,paste0(outcomevariable,"_hi")] <- {if(negative==FALSE){
      ifelse(output[,paste0(outcomevariable,"_predictedTE")]>=
               quantile(unique(output[,paste0(outcomevariable,"_predictedTE")]),c(1,2,3,4,5)/5)[4],1,0)
    }else{
      ifelse(output[,paste0(outcomevariable,"_predictedTE")]<=
               quantile(unique(output[,paste0(outcomevariable,"_predictedTE")]),c(1,2,3,4,5)/5)[1],1,0)
    }}
  }else{
    output <- cbind(hte_effect)%>%
      `colnames<-`(paste0(outcomevariable,"_predictedTE"))%>%
      as.data.frame
  }
  
  # put the results in the result list
  # NOTE: "<<-" means this result is forced to be saved in the global environment
  # and this step is set to make sure this function could be used for many different outcome variables
  # we do not need to rewrite the function when a new outcome variable comes in
  results[[paste0(outcomevariable,"_causalTree_matching_",length(results)+1)]] <<- output
  results[[paste0(outcomevariable,"_causalTree_matching_Tree_",length(results)+1)]] <<- opfit
  results[[paste0(outcomevariable,"_causalTree_matching_ttable_",length(results)+1)]] <<- ttable
}


# :::::::::::::::::::::::::::::::::::::::::
# feed grf estimates into causal tree leaves
# :::::::::::::::::::::::::::::::::::::::::
# NOTE: need to first run grf algorithms
hte_forest <- function(outcomevariable, 
                       # the name of outcome variabls we are interested in
                       minsize=20, 
                       # minimum number of treated observations, 
                       # control observations in a leaf
                       crossvalidation = 20, 
                       # number of cross-validations to do
                       data = edurose_mediation_20181126, 
                       # can be changed, and the defaul one defined here 
                       # is edurose_mediation_20181126, the education dataset we are
                       # working on
                       treatment_indicator = 'compcoll25', # treatment variable
                       ps_indicator = 'propsc_com25', # propensity scores 
                       ps_linear = 'propsc_com25lin',
                       covariates = c(linear_terms,ps_indicator),
                       negative = FALSE, 
                       # can be changed, specify the expected direction 
                       # of the treatment effects
                       drawplot = TRUE, 
                       # export the graph of tree structure if true
                       con.num=1,
                       # the number of control variables used in matching
                       no_indicater="",
                       legend.x = 0.08,
                       legend.y = 0.25,
                       gf=tau.forest,
                       ...){
  
  if(!exists("tau.forest")){
  stop("please first run grf algorithm and name is as tau.forest")
   }
  # delete all missing values which is required in causal tree model 
  # and use it as the train set in machine learning model
  trainset <-  data[!is.na(data[,outcomevariable]),]
  
  # set up the formula used for constructing causal tree
  # export the formula for causal tree model: Y~X
  if(nchar(ps_indicator)>0){
    covariates_ <- c(covariates,ps_indicator) # non-linear ps score
    formula <- as.formula(paste(outcomevariable," ~ ", 
                                paste(covariates,collapse = '+'), collapse= "+"))
    covariates <- c(covariates,ps_linear) # linear ps score
  }else{
    formula <- as.formula(paste(outcomevariable," ~ ", 
                                paste(covariates,collapse = '+'), collapse= "+"))
    covariates <- covariates
  }
  # set up propensity score
  if(nchar(ps_indicator)>0){
    # contruct tree
    tree <- causalTree(formula, 
                       # specify the model, outcome variable ~ covariates
                       data = trainset, # specify the dataset to be used
                       treatment = trainset[,treatment_indicator], 
                       # specify the treatment variable, must be 0/1 indicator
                       split.Rule = "CT", 
                       # specify split rule; for causal tree, use "CT" 
                       # NOTE: there are four different splitting rules, 
                       # they are different in the cross-validation criteria used
                       # to determine the tree structure
                       # 1 - TOT
                       # 2 - CT
                       # 3 - fit
                       # 4 - tstats
                       # 5 - totD
                       # 6 - ctD
                       # 7 - fitD
                       # 8 - tstatsD
                       cv.option = "CT", # specifify cross validation method 
                       # and there are four different methods -- tot, ct, fit, tstats
                       # for causal tree, use "CT"
                       split.Honest = T, cv.Honest = T, split.Bucket = F, 
                       
                       xval = crossvalidation, 
                       # number of cross-validations to do and the default number is 20
                       cp = 0, 
                       propensity = trainset[,ps_indicator],
                       # specify the propensity score; if is not specified, it will use sum(treatment) / nobs as the propensity score
                       minsize = minsize # minimum number of treated observations, control observations in a leaf
                       # the default minimum size is 20, according to Jennie and Yu Xie's paper (Estimating Heterogeneous Treatment Effects with Observational Data, 2012)
    )}else{
      tree <- causalTree(formula, 
                         data = trainset, 
                         treatment = trainset[,treatment_indicator], 
                         split.Rule = "CT", 
                         cv.option = "CT",
                         split.Honest = T, cv.Honest = F, split.Bucket = F, 
                         xval = crossvalidation, 
                         cp = 0, 
                         minsize = minsize 
      )
    }
  
  # prune this tree model to avoid the overfitting issues
  # get the complexity parameter (cp) to be trimmed--the least important splits
  opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
  # recursively snipping off the least important tree based on the complexity parameter (cp)
  opfit <- prune(tree, opcp)
  
  # matchin in leaves and return the predicted heterogeneous treatment effect
  prty <- as.party(opfit)
  opfit_tree <- as.Node(prty)
  
  hte_effect_setup <- list()
  # Matching Algorithms
  x_num <- 0
  opfit_tree$Do(function(node){
    x_num <<- x_num+1
    
    subgroup <- (trainset%>%rownames)%in%(node$data%>%rownames)
    print(subgroup%>%sum)
    x <- average_treatment_effect(gf,
                                  target.sample = "treated",
                                  subset = subgroup)
    
    z <- x[1]/x[2]
    p = exp(-0.717*z - 0.416*(z^2))
    
    hte_effect_setup[[x_num]] <<- cbind(x[1],p,x[2],
                                        nrow(node$data)/nrow(trainset)*100%>%round(.,1))
    
    # keep all related numbers in the environment of node
    node$predicted <- x[1]%>%as.numeric
    node$pvalue <- p
    node$standarderror <- x[2]%>%as.numeric
    node$samplesize <- nrow(node$data)
  })
  
  opfit_tree<<-opfit_tree
  hte_effect <- opfit_tree$Get("predicted")%>%as.numeric
  opfit$frame$yval <- hte_effect
  
  # create a new variable indicating the estimated treatment effect for each unit
  hte_effect <- opfit$frame$yval[opfit$where]
  
  # statistics
  ttable <<- unlist(hte_effect_setup)%>%matrix(.,ncol = 4,byrow = TRUE)%>%
    as.data.frame%>%
    `colnames<-`(c("Estimator","pvalue","se","SampleSize"))
  st <- rep("",length(ttable$pvalue))
  st[ttable$pvalue<0.05] <- "*"
  st[ttable$pvalue<0.01] <- "**"
  st[ttable$pvalue<0.001] <- "***"
  ttable$star <<- st
  # ttable$star <<- ifelse(ttable$pvalue<0.1,"*","")
  ttable$SampleSize <- round(ttable$SampleSize,1)
  
  # If makeing plots, the values from the original tree should be 
  # adjusted to the value generated from matching methods
  # adj_effect <- table(hte_effect)%>%as.data.table
  # opfit$frame$yval[match(adj_effect$N,opfit$frame$n)] <- as.numeric(adj_effect$hte_effect)
  
  
  # if drawplots is TRUE, make plots and export the plots
  if(drawplot==TRUE){
    # makeplots(opfit,gph,trainset,covariates,outcomevariable)
    makeplots(negative=negative, opfit.=opfit,gph=gph,trainset,
              covariates,outcomevariable,data.=data,ttable,no_indicater,legend.x,legend.y)
  }else{
    print(c('Drawplot = ', drawplot))
  }
  
  # In observational study:
  # export the estimation results by adding three new indicators:
  # *_hi variable is 0/1 variable indicating whether the treatment effect has excess the threshhold
  # *_predictedTE is the predicted treatment effect (but causal tree approach predicts the effects of a treatment theoretically)
  # *_rank is the ranking from highest to lowest relative to other units
  
  # merge all the variables into one single data frame
  # first put them in different columns, and add the ID: R0000100
  
  # In simulation study: just export the estimation results for treatment effects
  if(identical(trainset,edurose_mediation_20181126)){
    output <- cbind(out_hi=0,
                    out_predict=hte_effect,
                    R0000100 = data$R0000100) %>%
      # change it as a data table and rerank the predicted heterogeneous treatment effect
      as.data.table %>%
      # rerank
      .[order(hte_effect,decreasing = TRUE)] %>%
      # add the rankings relative to other observations
      .[,paste0(outcomevariable,"_rank"):=1:nrow(.)]
    # rename the column names for this new data frame
    colnames(output)[1:2] <- paste0(outcomevariable,c("_hi","_predictedTE"))
    # make sure it is a data frame
    output <- as.data.frame(output)
    # define the variable: *_hi which is a 0/1 variable indicating whether the treatment effect has excess the threshhold
    # And the threshhold is set to be the top 20% highest treatment effect
    output[,paste0(outcomevariable,"_hi")] <- {if(negative==FALSE){
      ifelse(output[,paste0(outcomevariable,"_predictedTE")]>=
               quantile(unique(output[,paste0(outcomevariable,"_predictedTE")]),c(1,2,3,4,5)/5)[4],1,0)
    }else{
      ifelse(output[,paste0(outcomevariable,"_predictedTE")]<=
               quantile(unique(output[,paste0(outcomevariable,"_predictedTE")]),c(1,2,3,4,5)/5)[1],1,0)
    }}
  }else{
    output <- cbind(hte_effect)%>%
      `colnames<-`(paste0(outcomevariable,"_predictedTE"))%>%
      as.data.frame
  }
  
  # put the results in the result list
  # NOTE: "<<-" means this result is forced to be saved in the global environment
  # and this step is set to make sure this function could be used for many different outcome variables
  # we do not need to rewrite the function when a new outcome variable comes in
  results[[paste0(outcomevariable,"_causalTree_matching_",length(results)+1)]] <<- output
  results[[paste0(outcomevariable,"_causalTree_matching_Tree_",length(results)+1)]] <<- opfit
  results[[paste0(outcomevariable,"_causalTree_matching_ttable_",length(results)+1)]] <<- ttable
}

