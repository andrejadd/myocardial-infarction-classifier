## Andrej Aderhold, 2016
##
##
##
## Define the place where all remaining scripts are located
##
## Note: The script for each classification method is sourced (loaded) only right befor
##       it is executed. This is to prevent premature loading in case there is a problem with
##       one of the scripts. Disable the corresponding method in the variable 'METHOD_LIST' below
##       to avoid loading of broken scripts. It it is broken, it is likely because an external package
##       can't be loaded. 
##



rm(list=ls())

SCRIPTS_DIR = "Scripts"

source(paste(SCRIPTS_DIR, "/sens_spec_error.R", sep=''))       ## sensitivity and specificity calculation
source(paste(SCRIPTS_DIR, "/quickhull_roc.R", sep=""))   ## convex hull of the ROC curve given fpr/sensitivity pairs, where fpr = 1 - specificity
source(paste(SCRIPTS_DIR, "/trapz.R", sep=""))           ## area under the ROC curve
source(paste(SCRIPTS_DIR, "/create_class_table.R", sep=""))           ## compiles Latex table with prediction scores for all methods.




## stop programm on error and recover in DEBUG mode
#options(error = recover)

##
## Specify classification methods to include. 
## Note: Methods are not executed in the order they appear here.
##
                                        
METHOD_LIST = c("GLM" , "LDA", "Decision-Tree", "KNN", "Lasso", "Uni-variate", "Random-Forest", "GP-ARD")
# METHOD_LIST = c("Decision-Tree")



##
## Setup test: 3 different setup of features included.
##
##
##   1 : Ta, Tref, SBP, ED, Ta/Tref, Ta/SBP
##   2 : Ta, Tref, SBP, ED 
##   3 : Ta, ED, Ta/Tref, Ta/SBP
##
DATA_SETUP_TEST = 3

##
## Number of times Lasso (and Elastic Net) are repeated with cross-validation runs
##  in order to determine a robust average of the regression coefficients.
##  The results between CV runs can slightly differ depending on how the random generator determines the training
##  and test sets for cross-validiation. This impacts the choice of the penalty parameter lambda.
##  In some instances, depending on the data, lambda can change significantly.
##  Repeating CV runs proved to be more stable than a single run but more data could make this step unnecessary 
##  and LASSO_BS_* could be set to 1 in the future.
##
LASSO_BS_PREDICTION = 20      ## A value of 20 is used in the paper, for prediction.
LASSO_BS_COEFF_ESTIMATE = 100 ## A value of 100 is used in the paper, for feature importance.


## Plotted results and .csv files go into this directory.
## Each method gets its own subdirectory.
RESULT_DIR = "Results"

## Enable or disable the importance plots
DO_IMPORTANCE_PLOTS = 1

## Turn on error rate calculation for different Random Forest parameter settings
## This takes time and is not really necessary
RF_ERROR_RATES = 0

##
## Read in data
##
## This is version 2 : 27 healthy and 11 MI patient records : n = 38 
##
all_data = read.table("Data/heart_data_v2.csv", sep=";", header=TRUE)


##
## Extract the features dependent on the setup
##

if (DATA_SETUP_TEST == 1)
    heart_data = data.frame(Type=as.factor(all_data$Type), Ta=all_data$Ta, Treq=all_data$Treq, SBP=as.numeric(all_data$SBP), EDV=all_data$EDV, "Ta/Treq"=all_data$Ta_Treq, "Ta/SBP"=all_data$Ta_SBP)

if (DATA_SETUP_TEST == 2)
    heart_data = data.frame(Type=as.factor(all_data$Type), Ta=all_data$Ta, Treq=all_data$Treq, SBP=all_data$SBP, EDV=all_data$EDV)

if (DATA_SETUP_TEST == 3)
    heart_data = data.frame(Type=as.factor(all_data$Type), Ta=all_data$Ta, EDV=all_data$EDV, "Ta/Treq"=all_data$Ta_Treq, "Ta/SBP"=all_data$Ta_SBP)



##
## Replace dot in name with devision symbol "/". The previous data.frame() converts "/" always into ".".
##
idtmp = which(names(heart_data) == "Ta.Treq")

if (length(idtmp) != 0)
    names(heart_data)[idtmp] = "Ta/Treq"

idtmp = which(names(heart_data) == "Ta.SBP")

if (length(idtmp) != 0)
    names(heart_data)[idtmp] = "Ta/SBP"

##
## Get some meta data from the data
##
variable_names = names(heart_data)[-1]           # get feature names
n_variables    = ncol(heart_data) - 1            # number of features/variables
n_samples      = nrow(heart_data)                # number of samples/observations
true_labels    = as.numeric(heart_data$Type) - 1 # get the true classes stored as 0 or 1 in the column with name "true"


## Convert variable string names to expression names. 
## This is for plotting subscripts and superscripts correctly.
## NOTE: Adapt these according to you needs. 
exp_var_names = c()   ## all expressions go into here, same size as 'variable_names'

## make some global substitutions 
varnamestmp = gsub("Ta", "T[a]", variable_names)   # put the a into subscript
varnamestmp = gsub("Treq", "T[req]", varnamestmp)  # put the "req" into subscript
varnamestmp = gsub("T\\[a\\]/SBP", "T[a]^norm", varnamestmp)   # this is for T[a]/SBP -> T[a]^norm, however, for some reason it will not detect the full string
varnamestmp = gsub("T\\[a\\]/T\\[req\\]", "C^s", varnamestmp) # this is for T[a]/T[req] -> C^s 

# previous calling convention:
#varnamestmp = gsub("T\\[a\\]/T\\[req\\]", "T[a]^req", varnamestmp) # this aims at T[a]/T[req] -> T[a]^req 


## now convert each of the above strings into an expression using parse()
for(i in 1:length(varnamestmp)) {
  exp_var_names = c(exp_var_names, parse(text=varnamestmp[i])) 
}



## Each method gets an id, which is used for plotting the right symbol.
## Just start with 1 or higher. 
method_id = 1


##
## The specificity and sensitivity values for each classification method are saved in this data frame.
## 
##  * 1st  column is a method id that determines the shape and color in the final ROC plot.
##  * 2nd column is the methods name.
##  * 3rd and 4th columns are the specificity and sensitivity respectively. 
##
## Note: The following notation is used to add new rows to this data frame:
##
##     roc_points[nrow(roc_points)+1, ] = list(a_numeric_value, "a_method_name", a_number, another_number)
##
## It's important _not_ to use rbind() or c(). rbind() will mess up the character column into a factor, and c() will coerce all elements to character strings.
## 
roc_points = data.frame(method_id=numeric(), method_name=character(), fpr=numeric(), sensitivity=numeric(), error=numeric(), stringsAsFactors=FALSE)


##
## Save the importance measures for each method and feature into here (if it is supported by a method)
##
##  importance_values is a g x p matrix that holds the feature importance values in each row (for g methods that provide importance measures) 
##  importance_names is an array of strings with the corresponding method names.
##  importance_descr is an array of strings with a description of the importance measure. 
##
## Note: Using a mixed data frame would be more concise, however, the data frames columns have to be allocated dynamically.
##       E.g.  setNames(replicate(5,numeric(0), simplify = F), letters[1:5])
##       does this. However, it has to be mixed with the method_name. Did not attempt this.
##
importance_values = matrix(nrow=0,ncol=n_variables)
importance_names = c()
importance_descr = c()   






## ===================================================
##
##
## START OF CLASSIFICATION METHODS.
##
##
##

cat("\nExecuting classifiers for data setup ", DATA_SETUP_TEST, " ..\n")



##
##
## Lasso sparse logistic regression 
##
##
method_name = "Lasso"
if(method_name %in% METHOD_LIST) {

    cat("Running ", method_name, "\n")
    source(paste(SCRIPTS_DIR, "/main_", method_name, ".R", sep=""))
    RESULT_DIR_METHOD = paste(RESULT_DIR, "/", method_name, sep="")

    ## create sub-directory in the results folder if not existing
    if(!file.exists(RESULT_DIR_METHOD))
        dir.create(RESULT_DIR_METHOD)

    PENALTY_MIXING = 1  # 1 sets it to Lasso, 0.5 to Elastic Net
    

    cat(" * Bootstrapping Lasso with ", LASSO_BS_PREDICTION, " trials")
    
    
    ## make prediction with penalized Lasso regression  
    cl = lasso_main_bs(heart_data, PENALTY_MIXING, LASSO_BS_PREDICTION)

    ## add to ROC plot
    roc_points[nrow(roc_points)+1, ] = list(method_id, method_name, (1 - Specificity(true_labels, cl)), Sensitivity(true_labels, cl), Error_rate(true_labels, cl))

    

    cat("\n * Determining regression coefficients with Lasso..\n")

    ##
    ## Get the feature importance measure of Lasso by extraction of the mean 
    ## regression coefficient over 100 bootstraps. 
    res = lasso_coeff_bs(heart_data, PENALTY_MIXING, LASSO_BS_COEFF_ESTIMATE)
  

    ## Extract importance measure:
    ## Absolute value of mean coefficients.
    ## Using this because it is not affected by the coefficient error.
    ## If we would have more data, we would prefer the z-statistic or
    ## the feature inclusion probability (see below).
    impvalues = as.vector(abs(res$mean_coeff))
    importance_descr  = c(importance_descr, "absolute mean coefficients")

    ## We could use the z-values as importance indicator.
    ## However this value has the problem that coefficients with high error
    ##  become unimportant. Some coefficients will have high values
    ##  in some instances and 0 values in others, depending on the cross-validation.
    ##  These high coefficient values should not disappear because of an high error. 
    ##  
    # impvalues = as.vector(abs(res$z_value))
    # importance_descr  = c(importance_descr, "coefficient z-value")
    
    ##
    ## Alternative measure: feature inclusion probability ("edge probability")
    ##  Similar to z-statistics: suffers from bad CV penalty optimization of Lasso.
    # impvalues = as.vector(res$edge_prob)
    # importance_descr  = c(importance_descr, "feature inclusion probability")
    
    
    ## save to importance data objects, evaluated later..
    importance_values = rbind(importance_values, impvalues)
    importance_names  = c(importance_names, method_name)
    
    
    
    ##
    ## Create a Lasso plot of the regularization path. 
    ##  Run the method several times (last argument) to show how the penalty parameter lambda
    ##  can drastically change in dependence of the randomness of the cross-validation.
    ##
    plot_lasso_coeffs(heart_data, RESULT_DIR_METHOD, DATA_SETUP_TEST, 3, PENALTY_MIXING)
  

}



##
##NOTE: THIS IS TURNED OFF, BECAUSE IT IS VERY SIMILAR TO LASSO (SEE THE PENALTY_MIXING PARAMETER)
##
## Elastic Net, sparse logistic regression 
##
##
method_name = "Elastic-Net"
if(method_name %in% METHOD_LIST) {

    cat("Running ", method_name, "\n")
    source(paste(SCRIPTS_DIR, "/main_", method_name, ".R", sep=""))
    RESULT_DIR_METHOD = paste(RESULT_DIR, "/", method_name, sep="")

    ## create sub-directory in the results folder if not existing
    if(!file.exists(RESULT_DIR_METHOD))
        dir.create(RESULT_DIR_METHOD)

    PENALTY_MIXING = 0.5  # 1 sets it to Lasso, 0.5 to Elastic Net
    

    cat(" * Bootstrapping Elastic Net with ", LASSO_BS_RUNS, " trials.\n")
    cat(" * Determine regression coefficients with Elastic Net..\n")

    res = lasso_coeff_bs(heart_data, PENALTY_MIXING, LASSO_BS_RUNS)

    ##
    ## Plot mean values and probabilities
    ##

    pdf(paste(RESULT_DIR_METHOD, "/elasticnet_mean_coeff_datatest", DATA_SETUP_TEST, ".pdf", sep=""), width=7, height=5, paper='special')
    barplot(res$mean_coeff, main="ElasticNet (CV) mean coefficients",
            xlab="variables", names.arg = variable_names)
    dev.off()

    pdf(paste(RESULT_DIR_METHOD, "/elasticnet_pos_mean_coeff_datatest", DATA_SETUP_TEST, ".pdf", sep=""), width=7, height=5, paper='special')
    barplot(res$pos_mean_coeff, main="ElasticNet (CV) absolute mean coefficients",
            xlab="variables", names.arg = variable_names)
    dev.off()

    pdf(paste(RESULT_DIR_METHOD, "/elasticnet_edge_prob_datatest", DATA_SETUP_TEST, ".pdf", sep=""), width=7, height=5, paper='special')
    barplot(res$prob_edge, main="ElasticNet (CV) edge probabilities",
            xlab="variables", names.arg = variable_names)
    dev.off()
}



##
##
## Uni-variate Logistic Regression with GLM
##
##
method_name = "Uni-variate"
if(method_name %in% METHOD_LIST) {
  
  cat("Running ", method_name, "\n")
  source(paste(SCRIPTS_DIR, "/main_GLM.R", sep=""))
  
  method_id = method_id + 1
  
  ##
  ## Extract each feature and the response and make a univariate prediction.
  ##
  for(i in 2:ncol(heart_data)) {
    
    tmp_data = heart_data[,c(1,i)]
    
    cl = main_glm(tmp_data)
    
    ## convert fraction to 0 or 1
    cl[cl > 0.5] = 1
    cl[cl != 1 ] = 0
    
    roc_points[nrow(roc_points)+1, ] = list(method_id, "Univariate LR", (1 - Specificity(true_labels, cl)), Sensitivity(true_labels, cl), Error_rate(true_labels, cl))
    
  }
}


##
##
## Multivariate Logistic Regression with GLM 
##
##
method_name = "GLM"
if(method_name %in% METHOD_LIST) {

    cat("Running ", method_name, "\n")
    source(paste(SCRIPTS_DIR, "/main_", method_name, ".R", sep=""))
    
    method_id = method_id + 1
    cl = main_glm(heart_data)
    
    ## use probability of response to sample class
    # cl = rbinom(length(cl), 1, cl)

    ## convert fraction to 0 or 1
    cl[cl > 0.5] = 1
    cl[cl != 1 ] = 0

    roc_points[nrow(roc_points)+1, ] = list(method_id, "Multivariate LR", (1 - Specificity(true_labels, cl)),Sensitivity(true_labels, cl), Error_rate(true_labels, cl))

    log_mod = glm(formula= Type ~ ., data=heart_data, family=binomial)
    ##browser()
    summary(log_mod)
}



##
##
## Principle Component Analysis (PCA)
##
## This method could be used for "feature generation", i.e. replacing the features by the principal components.
## I did not include this method here since there is only a small number of features of which a even smaller number
## is highly relevant. There seems to be no need for this method, although it might prove usefull in the future.
##
##
method_name = "PCA"
if(method_name %in% METHOD_LIST) {

    cat("Running ", method_name, "\n")
    source(paste(SCRIPTS_DIR, "/main_", method_name, ".R", sep=""))
    RESULT_DIR_METHOD = paste(RESULT_DIR, "/", method_name, sep="")

    ## create sub-directory in the results folder if not existing
    if(!file.exists(RESULT_DIR_METHOD))
        dir.create(RESULT_DIR_METHOD)

    method_id = method_id + 1

    ## res = main_glm(heart_data)

    

}




##
##
## KNN
##
##
method_name = "KNN"
if(method_name %in% METHOD_LIST) {

    cat("Running ", method_name, "\n")
    source(paste(SCRIPTS_DIR, "/main_", method_name, ".R", sep=""))

    RESULT_DIR_METHOD = paste(RESULT_DIR, "/", method_name, sep="")

    ## create sub-directory in the results folder if not existing
    if(!file.exists(RESULT_DIR_METHOD))
        dir.create(RESULT_DIR_METHOD)

    method_id = method_id + 1

    ## Find best parameter k with MSE_test using LOOCV
    ## Also does some plotting.
    k_param = main_KNN_getParams(heart_data, RESULT_DIR_METHOD, DATA_SETUP_TEST)

    ## execute KNN with the best estimated k
    cl = main_KNN(heart_data, k_param)

    ## add to ROC curve
    roc_points[nrow(roc_points)+1, ] = list(method_id, method_name, (1 - Specificity(true_labels, cl)),Sensitivity(true_labels, cl), Error_rate(true_labels, cl))
 

    
} # end KNN




##
##
## LDA 
##
##
method_name = "LDA"
if(method_name %in% METHOD_LIST) {
    
    cat("Running ", method_name, "\n")
    source(paste(SCRIPTS_DIR, "/main_", method_name, ".R", sep=""))
    
    method_id = method_id + 1
    cl = main_LDA(heart_data)

    roc_points[nrow(roc_points)+1, ] = list(method_id, method_name, (1 - Specificity(true_labels, cl)), Sensitivity(true_labels, cl), Error_rate(true_labels, cl))
    
}



##
##
## Decision Tree
##
##
method_name = "Decision-Tree"
if(method_name %in% METHOD_LIST) {

    
    cat("Running ", method_name, "\n")
    source(paste(SCRIPTS_DIR, "/main_", method_name, ".R", sep=""))
    RESULT_DIR_METHOD = paste(RESULT_DIR, "/", method_name, sep="")
    
    ## create sub-directory in the results folder if not existing
    if(!file.exists(RESULT_DIR_METHOD))
        dir.create(RESULT_DIR_METHOD)

    method_id = method_id + 1

    ## Try out a couple of parameters settings and plot sensitivities, specificities, and error rates. 
    ## Find best boosting trials/cost penalty combination that minimizes the misclassification error rate using LOOCV.
    ## Try these boosting trials:
    boost_trials = c(1,10,20,30)
    max_penalty  = 15

     
    ## Instead, every time we predict the class of a single test observation (given the remaining p-1 training observations),
    ##  We run a full LOOCV on these p-1 training observations with different parameter settings.
    res = main_tree(heart_data, boost_trials, max_penalty)

    cl = res$predictions

    text_out_file = paste(RESULT_DIR_METHOD, "/boosts_for_each_prediction_datatest", DATA_SETUP_TEST, ".csv", sep="")
    write.csv(file=text_out_file, x=res$boosts)

    text_out_file = paste(RESULT_DIR_METHOD, "/costs_for_each_prediction_datatest", DATA_SETUP_TEST, ".csv", sep="")
    write.csv(file=text_out_file, x=res$costs)
    

    ## Add to final result
    roc_points[nrow(roc_points)+1, ] = list(method_id, method_name, (1 - Specificity(true_labels, cl)), Sensitivity(true_labels, cl), Error_rate(true_labels, cl))

    
    
    ## Make the plot and return parameters with lowest error. 
    ## Note: We could use these parameters to make a prediction for a _new_ data point that is not 
    ##       in the dataset 'heart_data'.
    ##       However, we don't use it to make a predictions, because each test point we used for validation, 
    ##       is also part of the training set, thus our prediction is over-optimistic 
    params = main_tree_getBestParams(heart_data, boost_trials, max_penalty, DATA_SETUP_TEST, RESULT_DIR_METHOD)    
    
    ##
    ## Plot tree with best error setting
    ##
    
    fileout = paste(RESULT_DIR_METHOD, "/tree_boost", params$best_boost, "_error", params$best_cost_penalty, "_dataset", DATA_SETUP_TEST,".pdf", sep="" )
    pdf(fileout, width=5, height=5, paper='special' )
    
    tmod = C5.0(Type ~ ., data = heart_data, trials = params$best_boost, error = params$best_cost_penalty)
    plot(tmod, trial=0)
    dev.off()

    ## Get the importance as usage of observations classified by a feature
    ## metric = "usage" or "splits"
    importance = t(C5imp(tmod, metric="usage"))  ## usage gives values in percent
    
    ## The values are out of order, need to extract each given the original order
    impvalues = c()  ## importance values in new order go into here
    
    for (i in 1:length(variable_names)) {
        idx = which(colnames(importance) == variable_names[i])
        
        impvalues = c(impvalues, importance[idx])
    }

    
    ## save to importance data objects, evaluated later..
    importance_values = rbind(importance_values, impvalues)
    importance_names  = c(importance_names, method_name)
    importance_descr  = c(importance_descr, "usage in percent")

    
}


##
##
## Random Forest
##
##
method_name = "Random-Forest"
if(method_name %in% METHOD_LIST) {

    cat("Running ", method_name, "\n")
    source(paste(SCRIPTS_DIR, "/main_", method_name, ".R", sep=""))
    RESULT_DIR_METHOD = paste(RESULT_DIR, "/", method_name, sep="")

    ## create sub-directory in the results folder if not existing
    if(!file.exists(RESULT_DIR_METHOD))
        dir.create(RESULT_DIR_METHOD)

    method_id = method_id + 1

    
    ##
    ## Try out different parameters:
    ##   mtrys - number of sub-selected features for each split
    ##   n_trees - size of forest
    ##
    ## This will plot n_trees against error rate for different feature sub sets.
    ## NOTE: No point in using the parameters given a lowest error rate because
    ##       the results will change anyway for each new run.
    ##
    if(RF_ERROR_RATES)
        main_randforest_tryParams(heart_data, RESULT_DIR_METHOD, DATA_SETUP_TEST)
    
    ## run method: Tried this before, but will give always different results. 
    ## Reason: randomness of feature selection for splits.
    #res = main_randforest(heart_data, params$mtrys, params$n_trees, T)
    
    ## Give no n_trees or mtry parameter: use default values
    ##FIXME: if main_randforest takes NULL arguments, handle it properly inside the function.
    res = main_randforest(heart_data, NULL, NULL, T, RESULT_DIR_METHOD)
    
 
    ##
    ## Evaluate individual tree predictions by taking majority vote
    ##
    predicti = res$indiv_class
    
    cl_indiv = apply(predicti, 2, mean)
    cl_indiv[cl_indiv > 0.5] = 1
    cl_indiv[cl_indiv != 1]  = 0
    
    cat(" * majority vote over individual trees, FPR: ", 1 - Specificity(true_labels, cl_indiv), ", Sensitivity: ", Sensitivity(true_labels, cl_indiv), "\n")
    
    ##
    ## Evaluate aggregated class probability
    # cl = rbinom(length(res$class_prob), 1, res$class_prob)
    cl = res$aggr_prob
    cl[cl > 0.5] = 1
    cl[cl != 1]  = 0
    
    cat(" * aggregated probabilities, FPR: ", (1 - Specificity(true_labels, cl)), ", Sensitivity: ", Sensitivity(true_labels, cl), "\n")
    
    roc_points[nrow(roc_points)+1, ] = list(method_id, method_name, (1 - Specificity(true_labels, cl)), Sensitivity(true_labels, cl), Error_rate(true_labels, cl))

    ## importance values
    impvalues = res$importance    
    
    ## save to importance data objects, evaluated later..
    importance_values = rbind(importance_values, impvalues)
    importance_names  = c(importance_names, method_name)
    importance_descr  = c(importance_descr, "mean decrease in accuracy")

}
    


##
##
## Gaussian Process with Automatic relevance determination (ARD)
##  Uses: Logistic Likelihood and squared exponential (SE) kernel with ARD.
##
##  Note: These results come from a .csv file that has been created with the Matlab script 'gp_ard.m'.
##
method_name = "GP-ARD"
if(method_name %in% METHOD_LIST) {

    cat("Extracting ", method_name, " results. [run gp_ard.m to update/produce results]\n")
    RESULT_DIR_METHOD = paste(RESULT_DIR, "/", method_name, sep="")
    
    method_id = method_id + 1

    ## the gp_ard.m script saves the specificity and sensitivity values into this file
    gp_result_file = sprintf("%s/gp_covSEard_likLogistic_infEP_dataset%i.csv", RESULT_DIR_METHOD, DATA_SETUP_TEST)

    ## check if it exists, if not you need to execute gp_ard.m
    if(!file.exists(gp_result_file)) {
    
        print("GP ARD result file not found: \n RUN SCRIPT gp_ard.m \n", gp_result_file)
    
    } else {
        cat(" * OK.\n")
        
        ## read results and add to roc list
        cl = as.matrix(read.table(gp_result_file, sep=","))
        
        roc_points[nrow(roc_points)+1, ] = list(method_id, method_name, (1 - Specificity(true_labels, cl)), Sensitivity(true_labels, cl), Error_rate(true_labels, cl))

     
    }


    ##
    ## Extract importance measures of GP-ARD
    ##
    infile = sprintf("%s/gpard_lengthscales_dataset%i.csv", RESULT_DIR_METHOD, DATA_SETUP_TEST)

    ## check if it exists, if not you need to execute gp_ard.m
    if(!file.exists(infile)) {
    
        cat(" ! Could not read GP-ARD importance file: ", infile, "\n")
    
    } else {
        cat(" * Read importance file.\n")
        
        ## read results and add to roc list
        impvalues = as.matrix(read.table(infile, sep=",", header=T))

        ## overwrite feature names (can differ in GP-ARD), otherwise rbind() will fail if colnames do not match
        colnames(impvalues) = variable_names  
        
         ## save to importance data objects, evaluated later..
        importance_values = rbind(importance_values, impvalues)
        importance_names  = c(importance_names, method_name)
        importance_descr  = c(importance_descr, "inverted relative length scales")

    }

    
    ##
    ## This will load multiple predictions with GP-ARD. Each prediction is from a binomial sampling given 
    ##  the GP posterior probability. However, taking a majority vote over these samples will not improve
    ##  the prediction but its pretty much the same as with using a decision boundary of 0.5 for a single GP-ARD execution.
    ##
    ## Note: this might fail if the data was not previously generated with gp_ard.m
    GP_USE_BINOMIAL_SAMPLES = 0


    if (GP_USE_BINOMIAL_SAMPLES) {
      
      cl_mat = matrix(nrow=0, ncol=n_samples)
      for(i in 1:20) {
  
          cltmp = read.table(sprintf("%s/gp_covSEard_likLogistic_infEP_dataset%i_trial%i.csv", RESULT_DIR_METHOD, DATA_SETUP_TEST, i), header=F, sep=",")
  
          cl_mat = rbind(cl_mat, cltmp)
          
      }
  
      ## make majority vote, is equal to taking the mean and using a 0.5 decision boundary (as shown here)
      cl = apply(cl_mat, 2, mean)
      cl[cl > 0.5] = 1
      cl[cl != 1 ] = 0
      
      roc_points[nrow(roc_points)+1, ] = list(method_id, method_name, (1 - Specificity(true_labels, cl)), Sensitivity(true_labels, cl), Error_rate(true_labels, cl))
    }

    
    ##
    ## Compare likelihoods and inference methods.
    ##
    
    error_rates_mat = NULL
    
    ## the gp_ard.m script saves the specificity and sensitivity values into this file
    for(data_setup_id in 1:3) {

      error_rates = c()
      
      ## Logit
      gp_result_file1 = sprintf("%s/gp_covSEard_likLogistic_infEP_dataset%i.csv", RESULT_DIR_METHOD, data_setup_id)
      gp_result_file2 = sprintf("%s/gp_covSEard_likLogistic_infLaplace_dataset%i.csv", RESULT_DIR_METHOD, data_setup_id)
      gp_result_file3 = sprintf("%s/gp_covSEard_likLogistic_infVB_dataset%i.csv", RESULT_DIR_METHOD, data_setup_id)
      
      ## read results and add to roc list
      cl1 = as.matrix(read.table(gp_result_file1, sep=","))
      cl2 = as.matrix(read.table(gp_result_file2, sep=","))
      cl3 = as.matrix(read.table(gp_result_file3, sep=","))
      
      error_rates = c(error_rates, Error_rate(true_labels, cl1))
      error_rates = c(error_rates, Error_rate(true_labels, cl2))
      error_rates = c(error_rates, Error_rate(true_labels, cl3))
      
      ## Probit
      gp_result_file1 = sprintf("%s/gp_covSEard_likErf_infEP_dataset%i.csv", RESULT_DIR_METHOD, data_setup_id)
      gp_result_file2 = sprintf("%s/gp_covSEard_likErf_infLaplace_dataset%i.csv", RESULT_DIR_METHOD, data_setup_id)
      #gp_result_file3 = sprintf("%s/gp_covSEard_likErf_infVB_dataset%i.csv", RESULT_DIR_METHOD, DATA_SETUP_TEST)
      
      ## read results and add to roc list
      cl1 = as.matrix(read.table(gp_result_file1, sep=","))
      cl2 = as.matrix(read.table(gp_result_file2, sep=","))
      #cl3 = as.matrix(read.table(gp_result_file3, sep=","))
      
      error_rates = c(error_rates, Error_rate(true_labels, cl1))
      error_rates = c(error_rates, Error_rate(true_labels, cl2))
      
      if(is.null(error_rates_mat))
        error_rates_mat = matrix(nrow=0, ncol=length(error_rates))
      
      error_rates_mat = rbind(error_rates_mat, error_rates)
      
      rownames(error_rates_mat)[data_setup_id] = sprintf('Data %i', data_setup_id)

    }

    colnames(error_rates_mat) = c("EP-Logit", "Laplace-Logit", "VB-Logit", "EP-Probit", "Laplace-Probit") 
    error_rates_mat = t(error_rates_mat)
    
    ## do the simple barplot and put the legend outside the plot to the topright 
    pdf(paste(RESULT_DIR_METHOD, "/barplot_setups.pdf", sep=""), width=6, height=5, paper='special')
    par(mfrow=c(1,1), mar=c(5,5,4,8))
    barplot(error_rates_mat, main = sprintf("%s Setups", method_name), ylab = "error rate", beside=TRUE, 
            legend.text = TRUE,
            args.legend=list(x = "topright", 
                             bty = "n",
                             inset = c(-0.4,0)
                             )
            ) 
    
   
    dev.off()
    
}





cat("Classification Done.\n")


##
## Create Feature Importance plots for selected methods and a summary importance rank plot.
##
## Needs:
##   importance_values : matrix with measures for each feature (row) and method (column)
##   importance_names  : vector with method names matching the row entries in importance_values
##   importance_descr  : vector with description of the used importance measure used by each method 
##   variable_names    : vector with feature names matching the column entries in importance_values
##





if(DO_IMPORTANCE_PLOTS) {

    cat("Start feature importance plots.\n")


    ## store ranks here
    rank_mat = matrix(nrow=0, ncol=ncol(importance_values))

    
    ##
    ## Plot mean values and probabilities
    ##

    for ( i in 1:nrow(importance_values)) {

        method_name = importance_names[i]
        imp_values = importance_values[i,]
        description = importance_descr[i]

        RESULT_DIR_METHOD = paste(RESULT_DIR, "/", method_name, sep="")
        
        cat(" * Plotting importance measure for ", method_name, "\n")
        
        ## do the simple barplot 
        pdf(paste(RESULT_DIR_METHOD, "/feature_importance_barplot_datatest", DATA_SETUP_TEST, ".pdf", sep=""), width=6, height=5, paper='special')
        barplot(imp_values, main = sprintf("%s", method_name), ylab = description, names.arg = exp_var_names, cex.axis = 1.4)
        dev.off()

        ##
        ## Get rank,highest values have highest rank
        ##
        
        ## Special treatment for Decision-Tree? 
        ## features with a usage == 0, get the lowest rank 1, otherwise they would appear important
        #impvalues[impvalues_tmp == 0] = 1
        
        rank_mat = rbind(rank_mat, rank(imp_values))
        rownames(rank_mat)[nrow(rank_mat)] = method_name
        
        
        
    }


    ##
    ## Plot mean rank
    ##

    rank_mat = rank_mat - 1

    pdf(paste(RESULT_DIR, "/feature_importance_mean_ranks_datatest", DATA_SETUP_TEST, ".pdf", sep=""), width=6, height=5, paper='special')
  
    ## do the simple barplot 
    par(mfrow=c(1,1), mar=c(5,5,4,8))
    barplot(rank_mat, main = "Feature Importance Summary", ylab = "cumulative rank", beside=FALSE, 
            legend.text = TRUE,
            args.legend=list(x = "topright", 
                             bty = "n",
                             inset = c(-0.43,0)
            ),
            names.arg = exp_var_names
    ) 
    dev.off()
    

}




## 
## Calculate the convex hull, i.e. the ROC curve
##
##
hull = quickhull_roc( roc_points$fpr, roc_points$sensitivity)

fpr = hull$x
sensit = hull$y


cat("Start plotting final ROC ...\n")


##
## plot ROC curve and coordinates for each method
##
roc_out_file = paste(RESULT_DIR, "/roc_curve_datatest", DATA_SETUP_TEST, ".pdf", sep="")
pdf(roc_out_file, width=5, height=5, paper='special' )

plot(fpr, sensit, col = "black",type = "l", lwd = 1, lty = 2, xlim = c(0,1), ylim = c(0,1), ylab = "Sensitivity", xlab = "1 - Specificity")




##
## plot individual method points
##
for(idx in 1:nrow(roc_points)) {

    ## use method id to select color and point symbol
    m = roc_points$method_id[idx]

    ## get scores and add some jitter
    x = roc_points$fpr[idx] + rnorm(1,0,0.005)
    y = roc_points$sensitivity[idx] + rnorm(1,0,0.005)  

                                        #points(x,y, col = m, lty=2, lwd=2, pch = m, cex = 1.5)
    points(x,y, col = 1, lty=2, lwd=2, pch = m, cex = 1.5)
}


##
## Create title for the different data setups
##


## make title that includes all the feature names in a single row
AUROC_score =  trapz(fpr, sensit)

title(sprintf("ROC, AUROC: %.2g", AUROC_score))

##
## Get unique names and method ids for the legend
##
method_ids = unique(roc_points$method_id)
method_names = unique(roc_points$method_name)


legend(0.45, 0.6, legend=method_names,
       pch = method_ids, pt.cex=1.2, pt.lwd=2,
       col = 1)
       ## col = method_ids )

dev.off() ## save as pdf


##
## Write results to .CSV file
##

text_out_file = paste(RESULT_DIR, "/roc_points_results_datatest", DATA_SETUP_TEST, ".csv", sep="")
write.csv(file=text_out_file, x=roc_points)

cat(" * Prediction table written to ", text_out_file)
cat(" * ROC plot saved to ", roc_out_file, "\n")


##
## This script creates an overview of prediction accuracy for all methods for all
##  three data sets. It uses the above .CSV files roc_points_results_dataset*.csv
##
## Before running it will check if all files exists.
##
file1 = paste(RESULT_DIR, "/roc_points_results_datatest1.csv", sep="")
file2 = paste(RESULT_DIR, "/roc_points_results_datatest2.csv", sep="")
file3 = paste(RESULT_DIR, "/roc_points_results_datatest3.csv", sep="")

if(file.exists(file1) && file.exists(file2) && file.exists(file3) ) {
    fileout = create_class_table()
    cat(" * Latex table with method overview was written to: ", fileout, "\n")
} else {
    cat(" !  Latex overview table not create - Run this script for all three data sets.\n")
}



