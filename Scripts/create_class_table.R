##
##
## Function create_class_table() creates a Latex tabular with the classification scores for all methods
##  given the roc_point .csv files in the Results directory
##


    
create_class_table <- function() {

    ## where the roc_points .csv files are stored
    RESULT_DIR = "Results"

    latex_str = "\\begin{tabular}{l l | c | c | c | c | c | c | c | c } \n \\hline \n & "
    
    ## CREATE HEADER
    roc_points = read.table(sprintf("%s/roc_points_results_datatest1.csv", RESULT_DIR), sep=",", header=T)

    ## flag for including only one univariate
    UNIVARIATE_INCLUDED = FALSE
    
    for(i in 1:nrow(roc_points)) {

      if (roc_points$method_name[i] == "Univariate LR") {
        if(UNIVARIATE_INCLUDED) 
          next
      
        UNIVARIATE_INCLUDED = TRUE
      }
      
#      if ((roc_points$method_name[i] == "Univariate LR") && (!UNIVARIATE_INCLUDED)) {
        ## flag as True, next time it will skip the univariate
#        UNIVARIATE_INCLUDED = TRUE
#      }
        
      latex_str = sprintf("%s & %s", latex_str, roc_points$method_name[i])
    }
    
    latex_str = sprintf("%s \\\\ \n", latex_str)
    
    
    DATA_SETUP_TEST = 1
    latex_str = sprintf("%s \\hline \n %s", latex_str, compile_ROC_latex_table(DATA_SETUP_TEST, RESULT_DIR) )

    DATA_SETUP_TEST = 2
    latex_str = sprintf("%s \\hline \n %s", latex_str, compile_ROC_latex_table(DATA_SETUP_TEST, RESULT_DIR) )

    DATA_SETUP_TEST = 3
    latex_str = sprintf("%s \\hline \n %s", latex_str, compile_ROC_latex_table(DATA_SETUP_TEST, RESULT_DIR) )


    ## INSERT END
    latex_str = sprintf("%s\\hline \n\\end{tabular} \n", latex_str)

    fileout = sprintf("%s/class_table.tex", RESULT_DIR)
    fileConn<-file(fileout)
    writeLines(latex_str, fileConn)
    close(fileConn)

    return(fileout)


}



compile_ROC_latex_table <- function(DATA_SETUP_TEST, RESULT_DIR) {
    
    
    
    roc_points = read.table(sprintf("%s/roc_points_results_datatest%i.csv", RESULT_DIR, DATA_SETUP_TEST), sep=",", header=T)
    
    
    ##
    ## Find univariate with lowest error , this will be put into the 
    ##  into the table. The others will be ignored.
    ##
    
    min_error = 1000
    univariate_idx = -1
    
    for(i in 1:nrow(roc_points)) {
          
      if (roc_points$method_name[i] == "Univariate LR")  {
        ## check if error is smaller than before
        if (roc_points$error[i] < min_error) {
          min_error = roc_points$error[i]
          univariate_idx = i  ## remember index of the univariate
        }
      }
    }
    
    ##
    ## add Error
    ##
    
    la_str_tmp = sprintf("  & Error ")

    UNIVARIATE_INCLUDED = FALSE
    
    for(i in 1:nrow(roc_points)) {

        if (roc_points$method_name[i] == "Univariate LR")  {
          if(i != univariate_idx)
            next
          
          ## flag as True, next time it will skip the univariate
          ## UNIVARIATE_INCLUDED = TRUE
        }
         
        la_str_tmp = sprintf("%s & %.2g", la_str_tmp, roc_points$error[i])
        
    }
    la_str_tmp = sprintf("%s \\\\ \n", la_str_tmp)
    
    ##
    ## add Specificity
    ##
    
    ## flag for including only one univariate
    UNIVARIATE_INCLUDED = FALSE
    
    la_str_tmp = sprintf("%s $D_%i$ & Specificity ", la_str_tmp, DATA_SETUP_TEST)
    for(i in 1:nrow(roc_points)) {

        if (roc_points$method_name[i] == "Univariate LR"){
          if(i != univariate_idx)
            next
          
          #if(UNIVARIATE_INCLUDED)
          #  next
          
          ## flag as True, next time it will skip the univariate
          #UNIVARIATE_INCLUDED = TRUE
        }
      
        la_str_tmp = sprintf("%s & %.2g", la_str_tmp, (1 - roc_points$fpr[i]))
        
    }
    la_str_tmp = sprintf("%s \\\\ \n", la_str_tmp)
    
    ##
    ## add Sensitivity
    ##
    
    ## flag for including only one univariate
    UNIVARIATE_INCLUDED = FALSE
    
    la_str_tmp = sprintf("%s & Sensitivity ", la_str_tmp)
    for(i in 1:nrow(roc_points)) {

        if (roc_points$method_name[i] == "Univariate LR") {
          if(i != univariate_idx)
            next
          
          #if(UNIVARIATE_INCLUDED)
          #  next
          
          ## flag as True, next time it will skip the univariate
          #UNIVARIATE_INCLUDED = TRUE
        }
         
        la_str_tmp = sprintf("%s & %.2g", la_str_tmp, (roc_points$sensitivity[i]))
        
    }
    la_str_tmp = sprintf("%s \\\\ \n", la_str_tmp)
        

    return(la_str_tmp)
}
