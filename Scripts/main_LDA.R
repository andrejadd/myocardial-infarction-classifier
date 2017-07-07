
require(MASS)


##
## Linear Discriminant Analysis
##
main_LDA <- function(data2) {


    results = list()
      
    nr_samples = nrow(data2)    
    true.classes = data2[,1]

    
    ##
    ## Run LOOCV on LDA
    ##
    
    
    plda2preds <- c()
    for(i in 1:nr_samples){
        
       
        trainingset <- data2[-i,]

        train_cl = trainingset[,1] ## the class labels
        train = trainingset[,-1]   ## exclude class

        test  = data2[i,]
        test  = test[-1]  ## exclude class

        z = lda(train, train_cl)
        plda2preds = c(plda2preds, predict(z, test)$class)
        

    }

    ## save into results structure
    prediction = plda2preds - 1

    return(prediction)

}
