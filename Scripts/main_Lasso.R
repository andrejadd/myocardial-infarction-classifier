
require(glmnet)


##
## Repeat Lasso 'nr_bs' number of times and average over the predictions
##
lasso_main_bs <- function(data2, PENALTY_MIXING, nr_bs) {
    
    
    prediction_mat = matrix(0,ncol=nrow(data2), nrow=0)
    
    for (i in 1:nr_bs) {
        
        prediction = lasso_main(data2, PENALTY_MIXING)
        cat(".")
        prediction_mat = rbind(prediction_mat, prediction)
        
    }

    ## Use majority vote over bootstraps. This is equivalent to
    ## taking the mean and using a decision boundary of 0.5.

    mean_pred = apply(prediction_mat, 2, mean)

    ## Use decision boundary to select 0 or 1
    mean_pred[mean_pred > 0.5] = 1
    mean_pred[mean_pred != 1 ] = 0

    ## Alternative: sample from binomial given the probability
    # mean_pred = rbinom(length(mean_pred), 1, mean_pred)

    return(mean_pred)
}

##
## Make Lasso predictions with nested 10-fold cross-validation inside LOOCV
##
lasso_main <- function(data2, PENALTY_MIXING) {

    nr_samples = nrow(data2)

    predictions = c()

    results = list()
    
    ##
    ## Use LOOCV 
    ##

    plda2preds <- c()
    for(i in 1:nr_samples){

        trainingset <- data2[-i,]

        train_cl = as.vector(trainingset[,1])
        train = as.matrix(trainingset[,-1])

        test  = data2[i,]
        test  = as.matrix(test[-1])  ## exclude type

        ## run cross-validation on the training set
        cvres = cv.glmnet(train, train_cl, family="binomial", alpha = PENALTY_MIXING)

        ## make model with best lambda value
        cv.model = glmnet(train, train_cl, family="binomial", lambda=cvres$lambda.min, alpha= PENALTY_MIXING)

        ## predict response for test sample with Lasso model
        cl_prob <- predict(cv.model, newx = test, s = cv.model$lambda.min, type="response")
        
        ## append predicted value
        predictions = c(predictions, cl_prob)

    }

    return(predictions)
}



plot_lasso_coeffs <- function(indata, plot_dir, DATA_SETUP_TEST, run_trials, PENALTY_MIXING) {


    X = as.matrix(indata[,2:ncol(indata)])
    y = as.vector(indata[,1])

    ##
    ## Run cross-validiation several times to get (maybe) different estimates
    ##  for a lambda penalty. It can happen that the lambdas are very different, and so
    ##  so are the corresponding regression coefficients. 
    ##
    min_log_lambda = c()
    for (run_trial in 1:run_trials) {
      
      ## run cross-validation
      cvres = cv.glmnet(X, y, family="binomial", standardize=TRUE, alpha=PENALTY_MIXING)
      res = coef(cvres, s = "lambda.min")
    
      min_log_lambda = c(min_log_lambda, log(cvres$lambda.min))
    }
          
     
    
    ## For plotting: run full regularization path, plot and insert lambda mark
    lasso_res = glmnet(X, y, family="binomial")


    
    cat(" * Plotting regularization path given ", run_trials, " Lasso cross-validations.\n")
    
    pdf(paste(plot_dir, "/lasso_lambda_set_dataset",  DATA_SETUP_TEST, ".pdf", sep=""), width=7, height=5, paper='special')
    
    ## plot path
    plot(lasso_res, "lambda", lwd=2)
    title("Lasso regularization path")
    abline(v=min_log_lambda, lty=2, lwd=2)
    
    dev.off()

    return(res)
}



##
## Estimate the regression coefficients with bootstrapping.
##
##
lasso_coeff_bs <- function(indata, penalty_alpha_mixing, n_times) {


    X = as.matrix(indata[,2:ncol(indata)])
    y = as.vector(indata[,1])

    ## the coefficients are saved into here
    betas = matrix(ncol=ncol(X), nrow=0)
    
    for (n in 1:n_times) {
    
        ## run cross-validation
        cvres = cv.glmnet(X, y, family="binomial", standardize=TRUE, alpha=penalty_alpha_mixing, nfolds=5)
        res = coef(cvres, s = "lambda.min")
        #browser()
        ## remove intercept coefficient
        res = res[-1]

        ## append to beta matrix
        betas = rbind(betas, res)

    }

    ## get mean for each coefficient
    mean_coeff = apply(betas, 2, mean)
    
    ## get variance over bootstraps
    std_dev_betas = sqrt(apply(betas, 2, var))
    
    ## calculate z-value: estimate of coefficient / std.error of coefficient
    z_value = mean_coeff / (std_dev_betas / sqrt(n_times)) 
    
    ## calculate probability that edge is non-zero
    betas[betas != 0] = 1   # set edges not zero to 1
    edge_prob = apply(betas, 2, sum) / n_times # sum up and divide by number of trials
    
    return(list(mean_coeff = mean_coeff, z_value = z_value, edge_prob = edge_prob))
}

