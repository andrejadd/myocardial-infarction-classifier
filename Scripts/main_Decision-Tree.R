
require(C50)



main_tree <- function(input_data, boost_trials, max_penalty) {
  
  nr_samples = nrow(input_data)
  
  cat(' * Making prediction for each data point ')
  
  ##
  ## Predict the class for each data point in input_data by finding out the best 
  ## parameter given the remaining training set. The remaining training LOOCV with tree prediciton
  ##
  
  predictions = c()
  
  ## for the record, save the boost and cost for each prediction
  boosts = c()
  costs = c()

  
  ## Outer LOOCV that makes the actual predictions for each observation
  for(i in 1:nr_samples){
    
    cat('.')
    ## extract training set by taking out one sample
    train = data.frame(input_data[-i,])
    
    ## Inner LOOCV: use treetraining for another LOOCV to find the best parameter (without the use of the test observation)
    params = main_tree_getBestParams(train, boost_trials, max_penalty)
     
    
    ## save the parameters for later
    boosts = c(boosts, params$best_boost)
    costs  = c(costs, params$best_cost_penalty)
    
    ## transform cost penalty to matrix 
    ## set the error, i.e. the cost penalty
    if(is.null(params$best_cost_penalty))
      errorcosts = NULL
    else
      errorcosts <- matrix(c(0,0,params$best_cost_penalty,0), nrow = 2) 
    
    ## use the parameter to predict the class for the test observation
    treefit = C5.0(Type ~ . , data = train, trials = params$best_boost, costs = errorcosts)
    
    ## extract test sample (excludes type in first column)
    test  = data.frame(input_data[i,2:ncol(input_data)])
    
    ## predict class
    treepredict <- predict(treefit, test, type = "class")
    
    predictions <- c(predictions, treepredict)
    
    
  }
  cat('\n')
  
  return(list(predictions = predictions - 1, boosts = boosts, costs = costs))
    
}



##
##
## This funcion tries out different boost settings and runs main_tree_getParams() for each.
##   This is to get the the error rates for each boost so we can make a single plot, each boost gets one line.
##
main_tree_getBestParams <- function(input_data, boost_trials, max_penalty, DATA_SETUP_TEST=NULL, OUT_DIR=NULL ) {


    min_error_rate  = 1000      ## save best error here
    best_cost_penalty   = -1    ## save associated cost penalty here
    best_boost = -1             ## save associated boost here       

    ## The penalty (x-axis) vs. error rate (y-axis) for each boost goes into here
    ##  While x_penalties is a single vector (it is the same for all boosts), the y_error_rates is 
    ##  a matrix, each row corresponds to one boosting setting.
    y_error_rates = NULL
    x_penalties = c()
    
    for(boost in boost_trials) {
    
      ##  This will plot different views on the error and sensitivity and specificity given varying cost penalties.
      ##  Also returns the lowest cost_penalty and the associated misclassification error. 
      params_tmp = main_tree_getParams(input_data, boost, max_penalty, DATA_SETUP_TEST, OUT_DIR)
      
      ## Check if the error is smaller as what we have seen before
      if(params_tmp$error_rate < min_error_rate) {
        min_error_rate = params_tmp$error_rate       ## if smaller, save it along with..
        best_cost_penalty = params_tmp$cost_penalty  ## .. the associated penalty
        best_boost = boost                           ## and the boosting trials
      }
      
      ## save the penalty values, gets overwritten each time but doesn't matter since it stays the same
      x_penalties = params_tmp$penalties
      
      ## create matrix on the fly
      if(is.null(y_error_rates))
        y_error_rates = matrix(nrow=0, ncol=length(params_tmp$all_error_rates))

      ## append the error rates for this boost
      y_error_rates = rbind(y_error_rates, params_tmp$all_error_rates)
      
      ## free memory and manually invoke garbage collection
      rm(params_tmp)
      gc()
    }

    
    ##
    ## Do plotting if the following two function arguments are not NULL
    ##
    if ( !is.null(DATA_SETUP_TEST) && !is.null(OUT_DIR)) {
      
      ##
      ## plot cost penalty against the misclassification error rate
      ##
      pdf(paste(OUT_DIR, "/error_vs_cost_dataset", DATA_SETUP_TEST, ".pdf", sep=""), width=5, height=5, paper='special' )
      
      y_max_value = max(y_error_rates)  ## get largest y
      
      ## plot the error rates for each boost - loop through matrix y_error_rates.
      for(i in 1:nrow(y_error_rates)) {
        ## plot line
        plot(x_penalties, y_error_rates[i,],  type="o", lty=i, lwd=2, xlab="", ylab="", ylim=c(0, y_max_value)) 
        par(new=T)
        
        ## create name for label
        # legend_label = c(legend_label, sprintf('boost trials'))
      }
      mtext("cost penalty", 1, line=3, cex=1.5)
      mtext("error rate", 2, line = 3, cex=1.5)
      
      abline(v = best_cost_penalty, lty = 4, col="black")
      
      legend(0, y_max_value, lwd = 2,boost_trials, lty=1:nrow(y_error_rates), bg="white")
      
     
      title_str = sprintf("C5.0 Decision Tree with different boosting trials")
      title(title_str)
      
      dev.off()
      
    }    
    
    ## release memory
    rm(input_data)
    ## cat("Memory: ", memory.size(), "\n")
    return(list(best_cost_penalty=best_cost_penalty, best_boost=best_boost))

}



main_tree_getParams <- function(input_data, boost_trials, MAX_PENALTY, DATA_SETUP_TEST=NULL, OUT_DIR=NULL ) {

    ##
    ## Find out a cost penalty (error) that maximizes sensitivity/specificity scores using LOOCV.
    ##   Run main_tree() with different settings of the error, and extract the one that gives
    ##   the highest scores in terms of the sum of these scores.
    ##

    ## extract class labels as numeric values
    true_labels = as.numeric(input_data$Type) - 1
    
    roc_points = matrix(ncol=5, nrow=0)  ## columns: FPR, Sensitivity, penalty, Boosting, Flag (best)
    
    best_penalty = 0   # set with the best error, use for final run
    best_score = 0     # remember the best score for the last error

    error_rate_vec = c()
    
    for (penalty_cost in 0:MAX_PENALTY) {

        ## a penalty of 0 should be threated like a NULL (no penalty)
        penalty_cost_arg = penalty_cost
        if(penalty_cost == 0)
            penalty_cost_arg = NULL

        ## run the Tree method (uses LOOCV) with the error
        cl = main_tree_LOOCV(input_data, boost_trials, penalty_cost_arg)
        
        sens = Sensitivity(true_labels, cl)
        spec = Specificity(true_labels, cl)
        FPR  = 1 - spec
        
        roc_points = rbind(roc_points, c(FPR, sens, penalty_cost, boost_trials, 0))
        
        ## check if the score improved
        if( (sens + spec) > best_score) {
            best_score = sens + spec
            best_penalty = penalty_cost
        }

        ## Calculate misclassification error rate
        error_rate_vec = c(error_rate_vec, sum((cl - true_labels)^2) / length(true_labels))
    }

    ## flag ROC point that has highest (sensitivity + specificity)
    roc_points[best_penalty + 1 ,5] = 1
    
    ##
    ## Find index that has lowest misclassification error rate
    ##
    lowest_error_rate = min(error_rate_vec)
    best_penalty = which(error_rate_vec == lowest_error_rate)[1] - 1  ## subtract -1 to shift (1,..,) to (0,..,)
    
    penalties = 0:(length(error_rate_vec)-1)
    ## print out error penalty
    # cat(" * best cost penalty estimated (Sens/Spec) to be ", best_penalty, " using ", boost_trials, " boosting trial(s).\n")

    
    ##
    ## Do plotting if the following two function arguments are not NULL
    ##
    if ( !is.null(DATA_SETUP_TEST) && !is.null(OUT_DIR)) {
        
      ## 
      ## Calculate the convex hull, i.e. the ROC curve
      ##
      hull = quickhull_roc(roc_points[,1], roc_points[,2])
  
      fpr = hull$x
      sensit = hull$y
  
      ##
      ## plot ROC curve and coordinates for each method
      ##
      roc_out_file = paste(OUT_DIR, "/roc_curve_tree_boosting", boost_trials, "_dataset", DATA_SETUP_TEST, ".pdf", sep="")
      pdf(roc_out_file, width=5, height=5, paper='special' )
      
      plot(fpr, sensit, col = "black",type = "l", lwd = 1, lty = 2, xlim = c(0,1), ylim = c(0,1), ylab = "Sensitivity", xlab = "1 - Specificity")
  
  
      ##
      ## plot individual points, each one corresponds to a penalty setting
      ##
      for(idx in 1:nrow(roc_points)) {
  
          if(roc_points[idx,5] == 1)
            next
        
          ## get scores and add some jitter
          x = roc_points[idx,1] + rnorm(1,0,0.005)  # FPR
          y = roc_points[idx,2] + rnorm(1,0,0.005)  # sensitivity
          
          points(x,y, col = 1, lty=2, lwd=2, pch = 1, cex = 1.5)
      }
  
      ## draw the most left/upper point that indicates best score
      max_id = which(roc_points[,5] == 1)  ## this is the row with the flag set to 1
      
      par(new=T)
      plot( roc_points[max_id, 1], roc_points[max_id, 2], pch=4, lwd=3, cex=1.5, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
  
      ## add text to this point
      text_x = roc_points[max_id, 1]
      text_y = roc_points[max_id, 2] - 0.05
  
      ## make some hard coded adaptation, change if needed
      if(text_x > 0.8) {
          text_y = text_y + 0.06
      } else {
          text_x = text_x + 0.13
      }
      
      text(text_x, text_y, sprintf("penalty=%i", best_penalty))
  
      ## Create title for the different data setups
      var_names_str = paste(variable_names, collapse=', ')
  
      ## make title that includes all the feature names in a single row
      AUROC_score =  trapz(fpr, sensit)
      #title(sprintf("C5.0: Boosting trials: %i, AUROC: %.2g\n%s", boost_trials, AUROC_score, var_names_str))
      title(sprintf("C5.0: Boosting trials: %i, AUROC: %.2g", boost_trials, AUROC_score))
      
      dev.off()
    
          
      ## cat(" * lowest misclassification error is ", lowest_error_rate, " at penalty ", best_penalty, " using ", boost_trials, " boosting trial(s).\n")
        
      ##
      ## plot cost penalty against the misclassification error rate
      ##
      pdf(paste(OUT_DIR, "/error_vs_cost_boost", boost_trials, "_dataset", DATA_SETUP_TEST, ".pdf", sep=""), width=5, height=5, paper='special' )

      
      plot(penalties, error_rate_vec,  type="o", lty=1, lwd=2, xlab="penalty", ylab="error rate") 

      abline(v = best_penalty, lty = 4, col="black")

      title_str = sprintf("C5.0 Decision Tree, boosting: %i", boost_trials)
      title(title_str)

      dev.off()
    
    }
    
    ## free memory
    rm(input_data)
    
    ## return boosting trials and error
    return(list(boosting_trials = boost_trials, cost_penalty = best_penalty, error_rate = lowest_error_rate, all_error_rates = error_rate_vec, penalties = penalties ))
    
}



##
## Execute C5.0 Decision Tree prediction with LOOCV using BOOSTING_TRIALS and COST_PENALTY as parameters
##
main_tree_LOOCV <- function(input_data, BOOSTING_TRIALS, COST_PENALTY) {

    
    nr_samples = nrow(input_data)
        

    ##
    ## LOOCV with tree prediciton
    ##
    
    treepredictions <- c()

    ## set the error, i.e. the cost penalty
    if(is.null(COST_PENALTY))
        errorcosts = NULL
    else
        errorcosts <- matrix(c(0,0,COST_PENALTY,0), nrow = 2) 


    ## execute LOOCV iterations
    for(i in 1:nr_samples){

        ## extract training set by taking out one sample
        treetraining = data.frame(input_data[-i,])
        
        ## build tree
        treefit = C5.0(Type ~ . , data = treetraining, trials = BOOSTING_TRIALS, costs = errorcosts)

             
        ## extract test sample (excludes type in first column)
        treetesting  = data.frame(input_data[i,2:ncol(input_data)])

        ## predict class
        treepredict <- predict(treefit, treetesting, type = "class")
        
        treepredictions <- c(treepredictions, treepredict)

    }

    return(treepredictions - 1)
    
}





