

This folder contains R and Matlab scripts that assess the performance of several classification models that predict if a patient is under the risk of a myocardial infarction (MI).

The main R script is called "run_MI_classifiers.R". This script executes all classification methods, which are stored in the directory "Scripts". All methods are implemented in R, with the exception of the Gaussian Process method. This method is implemented in the Matlab script "gp_ard.m" also located in this directory. It has to be executed separately, and before the "run_MI_classifiers.R" script. The "run_MI_classifiers.R" scripts reads the result from the this Matlab script and creates all the plot outputs.

Usage:
---------

1. Execute "gp_ard.m" with Matlab.

2. Execute "run_MI_classifiers.R", e.g. in a R shell type 
   > source("run_MI_classifiers.R") 

3. Compile the report file that includes plot output from the above scripts:
   $ pdflatex MI_Class_Report.tex

Note 1: Store your data file in the directory "Data/". The file should be in the format of a tabular-separated text file. 
   The first row in this file is the header with the names of the outcome (Type) and the names of the  parameter/features that will be used to predict the outcome.
   Each row represents one sample, i.e. the parameters of one patient. 
   The outcome variable in the first column can be either '0' (healthy patient), or '1' (patient with MI).
   The remaining columns should be numbers and correspond to the measured parameters. See the data file "Data/heart_ver2.csv" as an example.

Note 2: Edit "gp_ard.m" and "run_MI_classifiers.R" so that the correct file is read.

Note 3: If you wish to exclude a classification method, edit the method list variable 'METHOD_LIST' in the file "run_MI_classifiers.R".  



Results:
----------

The script "run_MI_classifiers.R" will write all results to the directory "Results/" with sub folders for each method. 
The main result is a summary of all classification methods illustrated as a ROC curve in the file "Results/roc_curve_datatest[1|2|3].pdf". 
In the current version there are 3 data tests that differ in the parameters they use for prediction:

"Results/roc_curve_datatest1.pdf"  (if 'DATA_SETUP_TEST = 1' in "run_MI_classifiers.R")
				   Includes the four basic features (Ta, Tref, SBP, ED) and the two ratios Ta/Tref and Ta/SBP.

"Results/roc_curve_datatest2.pdf"  (if 'DATA_SETUP_TEST = 2' in "run_MI_classifiers.R")
				   Includes only the four basic features (Ta, Tref, SBP, ED).

"Results/roc_curve_datatest3.pdf"  (if 'DATA_SETUP_TEST = 3' in "run_MI_classifiers.R")
				   Includes two basic features (Ta,ED), and the two ratios Ta/Tref and Ta/SBP. 
				   Tref and SBP have been eliminated from the feature set because they are already in the ratios.

Note: You need to set the variable 'DATA_SETUP_TEST', and run the script "run_MI_classifiers.R" separately for each of the tree settings.

In addition, several plots are created for individual classification methods, such as for KNN, Decision Trees, and the Lasso. 
These plots are compiled into the Latex file MI_Class_Report.tex. Use 'pdflatex MI_Class_Report.tex' to create the corresponding pdf file.

A text output of all the sensitivity and specificity values is saved into the file "Results/roc_points_results.csv". 


Issues / Bugs
---------------

See the file Known_Issues.org for TODO points and known warnings and errors. Some of them can be fixed, others can be ignored.

