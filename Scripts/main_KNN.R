
require(class)

##
## K-nearest neightbours classification
##


##
## Estimates the k = (1,...,n-1) with the lowest misclassification error rate using LOOCV.
##   Returns 'k' with lowest MSE.
##
main_KNN_getParams <- function(heart_data, OUT_DIR, DATA_SETUP_TEST) {
                              


    results = list()


    ##
    ## Scale data : Important for KNN whenever variables ranges differ greatly.
    ##              KNN uses a distance measure and needs variables in similar intervals.
    ##
    
    data.sc = data.frame(scale(heart_data[,-1]))
    data.sc = cbind(Type=heart_data$Type, data.sc)
      
    nr_samples = nrow(data.sc)

    ## convert factors to numeric [0,1]
    true_labels = as.numeric(heart_data$Type) - 1
        
    ##
    ## Run KNN with LOOCV and varying parameters k = (1,...,n-1)
    ## Calculate MSE_test for each iteration of LOOCV, take average at the end
    ## and select k with lowest MSE.
    ##  Note: MSE for binary classification is called misclassification error rate: (FN + FP) / n = 1/n sum_i^n(y_i - f^_i)^2 , where f^_i is the estimated class label for observation i, y_i is the true label and n is the number of total observations.
    ##

    ## try all the k = (1,...,n-1) where n is the number of observations
    max_k = nr_samples - 1

    ## Fill in the predicted class into the following matrix
    ## For each sample in the LOOCV one column, for each parameter 'k' of KNN one row
    predict_mat = matrix(nrow=max_k, ncol=nr_samples)

    ## misclassification error rates
    error_vec = c()
    sens_vec  = c()
    spec_vec  = c()
    
    datatmp = data.sc
    
    for(k in 1:max_k) {
        
        ## Do LOOCV: take a test sample x_i and a training set (x_1, ..., x_i-1, x_i+1, ..., x_n).
        ## Make prediction with test sample x_i and calculate test MSE_i (BS a.k.a misclassification error rate)
        ##
        for(i in 1:nr_samples){
            
            train    = datatmp[-i, -1]
            train_cl = datatmp[-i, 1]
            
            test     = datatmp[i,-1]
            
            predict.model <- knn(train = train, test = test, cl = train_cl, k = k, prob = TRUE, use.all = TRUE)
            
            ## fill in the predicted class
            predict_mat[k,i] = unclass(predict.model) - 1

            
        }

        ## calculate error/sensitivity/specificity and append to vector        
        error_vec = c(error_vec, Error_rate(true_labels, predict_mat[k,]))
        sens_vec = c(sens_vec, Sensitivity(true_labels, predict_mat[k,]))
        spec_vec = c(spec_vec, Specificity(true_labels, predict_mat[k,]))
        
    }
    

    ## find index that has lowest MSE, is equal to 'k'
    best_k_error = which(error_vec == min(error_vec))[1]

    ##
    ## plot parameter k against error rate
    ##
    cat(" * Plotting k against error/sensitivity/specificity\n")
    pdf(paste(OUT_DIR, "/knn_error_k_dataset", DATA_SETUP_TEST, ".pdf", sep=""), width=5, height=5, paper='special' )
    
    
    plot(1 - error_vec, type="l", lty=1, lwd=2, xlim=c(1,length(error_vec)), ylim=c(0,1), xlab="k", ylab="")
    par(new=TRUE)
    plot(sens_vec, type="l", lty=2, lwd=2, xlim=c(1,length(error_vec)), ylim=c(0,1), xlab="", ylab="")
    par(new=TRUE)
    plot(spec_vec, type="l", lty=3, lwd=2, xlim=c(1,length(error_vec)), ylim=c(0,1), xlab="", ylab="")
    
    abline(v = best_k_error, lty = 4, col="black")
    title("K-Nearest Neighbours")
    
    legend(x="bottomright", c("1 - error rate", "sensitivity", "specificity"), lty=c(1,2,3), lwd=2, bg="white")
    par(new=FALSE)
    dev.off()
    

    

    ##
    ## Calculate the sensitivities and specificities for each parameter k
    ##  and plot th ese against k
    ##

    EXTRA_CALC = 0
    if(EXTRA_CALC) {
    sens_max = 0
    spec_max = 0

    sens_vec = c()
    spec_vec = c()

    sensspec_max = 0
    
    for(k in 1:nrow(predict_mat)) {

        senstmp = Sensitivity(true_labels, predict_mat[k,])
        spectmp = Specificity(true_labels, predict_mat[k,])

        sens_vec = c(sens_vec, senstmp)
        spec_vec = c(spec_vec, spectmp)

        if ( (senstmp + spectmp) > sensspec_max) {
            sensspec_max = senstmp + spectmp
            best_k = k
            sens_max = senstmp
            spec_max = spectmp
        }

    }

    ##
    ## Plot parameter k = (1,...,n-1) against sensitivity and specificity
    ##
    cat(" * Plotting k against (sensitivity,specificity).\n")
    pdf(paste(OUT_DIR, "/knn_sens_spec_k_dataset", DATA_SETUP_TEST, ".pdf", sep=""), width=5, height=5, paper='special' )

    plot(sens_vec, type="l", lty=2, lwd=2, ylim=c(0,1), xlim=c(1,nrow(predict_mat)), xlab="", ylab="")
    par(new=T)
    plot(spec_vec, type="l", lty=3, lwd=2, ylim=c(0,1), xlim=c(1,nrow(predict_mat)), xlab="k", ylab="Sensitivity and Specificity")

    abline(v = best_k, lty = 4, col="black")
    title("K-Nearest Neighbours")
    text(best_k + 3, 0, "selected k")
    dev.off()
    }
    
    return(best_k_error)
    
}


main_KNN <- function(heart_data, k_param) {
                              

    results = list()

    ##
    ## Prepare data for following executions
    ##
    data.sc = data.frame(scale(heart_data[,-1]))
    data.sc = cbind(Type=heart_data$Type, data.sc)
    datatmp = data.sc      

    nr_samples = nrow(data.sc)

    ## convert factors to numeric [0,1]
    true_labels = heart_data$Type - 1


    ## predicted class labels go into here
    predicted = c()

    
    ##
    ## Do LOOCV: take a test sample x_i and a training set (x_1, ..., x_i-1, x_i+1, ..., x_n).
    ##
    for(i in 1:nr_samples) {
            
        train    = datatmp[-i, -1]
        train_cl = datatmp[-i, 1]
        
        test     = datatmp[i,-1]
        
        predict.model <- knn(train = train, test = test, cl = train_cl, k = k_param, prob = FALSE, use.all = TRUE)
        
        ## fill in the predicted class
        predicted = c(predicted, unclass(predict.model) - 1)
        
    }

    return(predicted)
    
}



    
    
    