
require(randomForest)


##
## Try out different parameter setups for 
##   mtrys - number of sub-selected features for each split
##   n_trees - size of forest
##
## This will plot n_trees against error rate for different feature sub sets.
##
##
main_randforest_tryParams <- function(input_data, OUT_DIR, DATA_SETUP_TEST) {
  
  nr_features = ncol(input_data) - 1

  true_labels = unclass(input_data$Type) - 1
  
  cat(" * Calculating error rate for different forest size and feature subsets: ")
  start_trees = 1
  max_trees = 500
  
  error_rates_mat = NULL
  tree_seq = seq(start_trees, max_trees, 10)

  # help to find the parameters with the lowest error
  best_n_trees = -1
  best_mtrys   = -1
  lowest_error = 10000
  

  ##
  ## Loop over all the forest sizes.
  ##  For each forest size try out different feature sub-selection sizes.
  ##
  for(trees in tree_seq) {
    cat(".")
    error_rates = c()
  
    for(m in 1:nr_features) {
    
      res = main_randforest(input_data, m, trees, F, OUT_DIR)
    
      ## use decision boundary
      cl = res$aggr_prob
      cl[cl > 0.5] = 1
      cl[cl != 1]  = 0

      ## Calculate error rate
      error_tmp = Error_rate(true_labels, cl)
      
      ## .. and save
      error_rates = c(error_rates, error_tmp)
      
      ## remember the parameters if the error was very low
      if(error_tmp < lowest_error) {
        lowest_error = error_tmp
        best_n_trees = trees
        best_mtrys = m
      }
      
    }

    if(is.null(error_rates_mat))
      error_rates_mat = matrix(0, nrow=0, ncol=nr_features)
      
    error_rates_mat = rbind(error_rates_mat, error_rates)
  }
  
  pdf(paste(OUT_DIR, "/mtrys_vs_error_dataset", DATA_SETUP_TEST, ".pdf", sep=""), width=6, height=5, paper='special')

  ## put feature size for subselection as sting into here - becomes the legend
  m_trys_legend = c()

  ## plot the errors for each feature sub-selection size  
  for(i in 1:nr_features) {
    ## create label for the legend
    m_trys_legend = c(m_trys_legend, sprintf("size %i", i))
    
    ## plot line 
    plot(tree_seq, error_rates_mat[,i], type="o", pch=i, xlab="", ylab="", ylim=c(min(error_rates_mat),max(error_rates_mat)))
    par(new=T)
  }
  
  mtext("number of trees", 1, line=3, cex=1.5)
  mtext("error rate", 2, line = 3, cex=1.5)

  legend(x="topright", m_trys_legend, pch=1:nr_features, bg="white")
  
  title_str = sprintf("Random Forests with varying\nfeature sub-selection size")
  title(title_str)
  
  abline(v = best_n_trees, lty = 4, col="black")
  dev.off()
  
  cat("\n * best forest size: ", best_n_trees, ", best feature subset size: ", best_mtrys, ", error rate is ", lowest_error, "\n")
  
  return(list(mtrys = best_mtrys, n_trees = best_n_trees, lowest_error=lowest_error))

}



main_randforest <- function(input_data, mtrys, n_trees, IMPORTANCE, OUT_DIR) {

    results = list()

    nr_samples = nrow(input_data)
    
    ##
    ## Random forest with LOOCV
    ##

    ## for the individual predition of each tree
    predictions = NULL
    
    ## for the aggregated prediction
    class_prob = c()
    
    for(i in 1:nr_samples){

        train = input_data[-i,-1]   # training set withough first column (class)
        train_cl = input_data[-i,1] # only classes of training set

        test = input_data[i,-1]     # test sample without class

        ## randomForest detects classification whenever y is a vector of factors
        ##  ntree : size of the forest, i.e. number of trees
        ##  mtry  : number of features to extract for a split operation (feature sub-selection)
        ## 
        
        if (!is.null(mtrys) && !is.null(n_trees))        
          fitrf <- randomForest(x = train, y = train_cl, mtrys=mtrys, ntree = n_trees, importance=F)
        else
          fitrf <- randomForest(x = train, y = train_cl, importance=F)
        
        #
        # This plots the nr_trees vs. errors using the plot.randomForest() function
        #
        #pdf(paste(OUT_DIR, "/rf_plot_", DATA_SETUP_TEST, ".pdf", sep=""), width=6, height=5, paper='special')
        #plot(fitrf)
        #dev.off()
        
        ## make prediction
        predictrf <- predict(fitrf, test, type = "prob", predict.all = TRUE)

       
        ## extract the individual predictions of each tree
        if(is.null(predictions))
          predictions = matrix(nrow = length(predictrf$individual), ncol = 0)
        
        
        predictions = cbind(predictions, as.vector(predictrf$individual))
        
        ## extract the aggregated probability
        class_prob = c(class_prob, 1 - unclass(predictrf$aggregate)[1])
    }

    
    ## Gets prediction rates for the randomForest 
    ## and save for the return
    results$indiv_class <- matrix(as.numeric(predictions), ncol = ncol(predictions))
    results$aggr_prob = class_prob
    
    ##
    ## Estimate the importance. Run again but now with full data set and importance=TRUE
    ##
    
    if(IMPORTANCE) {
    
      fitrf = randomForest(x = input_data[,-1], y = input_data$Type, importance=T)
      results$importance = fitrf$importance[,3]   ## 3rd column is "MeanDecreaseAccuracy"
    
    } else {
      results$importance = NULL
    }
    
    return(results)

    }
