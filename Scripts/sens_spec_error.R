
## Sensitivity: 
##
## If a person has a disease, how often will the test be positive (true positive rate)? 
##
## Put another way, if the test is highly sensitive and the test result is negative you can be nearly certain that they don’t have disease. 
## A Sensitive test helps rule out disease (when the result is negative). Sensitivity rule out or "Snout"
##
## Sensitivity= true positives/(true positive + false negative)
##
Sensitivity <- function(true, test) {

    #if isnull(test) {
    
    ## number of true positives
    TP = sum(true & test)

    ## number of false negatives
    FN = sum(true & !test)

    return(TP / (TP + FN))
}


## Specificity: 
##
## If a person does not have the disease how often will the test be negative (true negative rate)?
##
## In other terms, if the test result for a highly specific test is positive you can be nearly certain that they actually have the disease.
## A very specific test rules in disease with a high degree of confidence Specificity rule in or "Spin".
##
## Specificity=true negatives/(true negative + false positives)
##
Specificity <- function(true, test) {

    ## number of true negatives
    TN = sum(!true & !test)

    ## number of false positives
    FP = sum(!true & test)

    return(TN / (TN + FP))
}

##
## Misclassification error rate, i.e. rate of false classifications
##  This score is equivalent to the 2 class Brier Score and MSE with 0 and 1 response.
##
Error_rate <- function(true, test) {

    ## number of false positives
    FP = sum(!true & test)

    ## number of false negatives
    FN = sum(true & !test)

    error = (FN + FP) / length(true)
    
    return(error)

}
