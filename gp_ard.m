% Andrej Aderhold, 2016, GP-ARD on MI data.
%
% The main script is gp_ard(). It will run GP-ARD and extract the length
% scales, and plot the posterior probabilities in a contour plot.
% Furthermore it will make the class predictions for the standard GP setup
% and for other setups.
%
% All results are written in the subdirectory specified by 'OUT_DIR'.
%


%
% Main function, runs several routines.
%
function gp_ard()

    clear all, close all
    
    % load the GPstuff package
    currpath = pwd();
    cd gpml-matlab;
    startup;
    cd(currpath);
    
    % If set to 1, it will run different likelihoods and inference methods on the data 
    TEST_MORE_GP_SETUP = 1;

    % PREDICTION_METHOD can be 'DECISION_BOUNDARY' (DEFAULT) or 'BINOM_SAMPLING' 
    PREDICTION_METHOD = 'DECISION_BOUNDARY';

    % scores and plots are saved into this directory 
    OUT_DIR = 'Results/GP-ARD';

    
    %
    % Load the data file. Note that this file has "noNamesColumn" as
    % postfix. The Matlab importdata() function has a problem with strings
    % in the matrix, thus a matrix without the names was created.
    datafile = 'Data/heart_data_v2_noNamesColumn.csv';
    
   
    % Load the different datasets, each uses the same data file but has
    % different features in it.
    [heart_data1, vnames1] = load_data(datafile, 1);
    [heart_data2, vnames2] = load_data(datafile, 2);
    [heart_data3, vnames3] = load_data(datafile, 3);
    
    %
    % Run GP-ARD and extract the length scales. Make posterior probability
    % contout plot for the most important features.
    %
    PLOT_LSCALE_BARPLOT = 0;    % Let R do the plotting.
    gpard_lscales_probmesh(heart_data1, vnames1, OUT_DIR, PLOT_LSCALE_BARPLOT, @likLogistic, 1)
    gpard_lscales_probmesh(heart_data2, vnames2, OUT_DIR, PLOT_LSCALE_BARPLOT, @likLogistic, 2)
    gpard_lscales_probmesh(heart_data3, vnames3, OUT_DIR, PLOT_LSCALE_BARPLOT, @likLogistic, 3)
    
    %
    % For final evaluation
    %
    gp_core(heart_data1, @covSEard, @likLogistic, @infEP, sprintf('%s/gp_covSEard_likLogistic_infEP_dataset1.csv', OUT_DIR), PREDICTION_METHOD); 
    gp_core(heart_data2, @covSEard, @likLogistic, @infEP, sprintf('%s/gp_covSEard_likLogistic_infEP_dataset2.csv', OUT_DIR), PREDICTION_METHOD);
    gp_core(heart_data3, @covSEard, @likLogistic, @infEP, sprintf('%s/gp_covSEard_likLogistic_infEP_dataset3.csv', OUT_DIR), PREDICTION_METHOD);

    
    %
    % The following runs are for comparing different likelihoods (Probit
    % a.k.a. Erf and Logit) and inference methods
    %
    if TEST_MORE_GP_SETUP

        gp_core(heart_data1, @covSEard, @likLogistic, @infLaplace, sprintf('%s/gp_covSEard_likLogistic_infLaplace_dataset1.csv', OUT_DIR), PREDICTION_METHOD); 
        gp_core(heart_data2, @covSEard, @likLogistic, @infLaplace, sprintf('%s/gp_covSEard_likLogistic_infLaplace_dataset2.csv', OUT_DIR), PREDICTION_METHOD);
        gp_core(heart_data3, @covSEard, @likLogistic, @infLaplace, sprintf('%s/gp_covSEard_likLogistic_infLaplace_dataset3.csv', OUT_DIR), PREDICTION_METHOD);

        gp_core(heart_data1, @covSEard, @likLogistic, @infVB, sprintf('%s/gp_covSEard_likLogistic_infVB_dataset1.csv', OUT_DIR), PREDICTION_METHOD); 
        gp_core(heart_data2, @covSEard, @likLogistic, @infVB, sprintf('%s/gp_covSEard_likLogistic_infVB_dataset2.csv', OUT_DIR), PREDICTION_METHOD);
        gp_core(heart_data3, @covSEard, @likLogistic, @infVB, sprintf('%s/gp_covSEard_likLogistic_infVB_dataset3.csv', OUT_DIR), PREDICTION_METHOD);

        gp_core(heart_data1, @covSEard, @likErf, @infEP, sprintf('%s/gp_covSEard_likErf_infEP_dataset1.csv', OUT_DIR), PREDICTION_METHOD); 
        gp_core(heart_data2, @covSEard, @likErf, @infEP, sprintf('%s/gp_covSEard_likErf_infEP_dataset2.csv', OUT_DIR), PREDICTION_METHOD);
        gp_core(heart_data3, @covSEard, @likErf, @infEP, sprintf('%s/gp_covSEard_likErf_infEP_dataset3.csv', OUT_DIR), PREDICTION_METHOD);

        gp_core(heart_data1, @covSEard, @likErf, @infLaplace, sprintf('%s/gp_covSEard_likErf_infLaplace_dataset1.csv', OUT_DIR), PREDICTION_METHOD); 
        gp_core(heart_data2, @covSEard, @likErf, @infLaplace, sprintf('%s/gp_covSEard_likErf_infLaplace_dataset2.csv', OUT_DIR), PREDICTION_METHOD);
        gp_core(heart_data3, @covSEard, @likErf, @infLaplace, sprintf('%s/gp_covSEard_likErf_infLaplace_dataset3.csv', OUT_DIR), PREDICTION_METHOD);
    
    end
    
end




%
% Main GP-ARD code that uses LOOCV to predict class labels (out of sample).
%
function[spec, sens, predictions] = gp_core(heart_data, covfunc, likfunc, inffunc, outfile_csv, PREDICTION_METHOD)


  
    % for Sensitivity/Specificity calculation
    true_classes = heart_data(:,1);

    % number of samples and variables (including class)

    nr_samples = size(heart_data, 1);
    nr_vars    = size(heart_data, 2) - 1; % excluding response

    fprintf('number variables: %i, number samples: %i\n', nr_vars, nr_samples);

    % the predicted values go into here
    class_probs = [];


    %
    % LOOCV loop over the data set
    %
    for i = 1:nr_samples

        % remove sample i from data matrix
        trainx = heart_data(setdiff(1:size(heart_data,1), i), :);

        % get test data
        testx = heart_data(i, 2:end);
       
        % extract response and change labels
        y = trainx(:,1);   % class response
        y(y == 0) = -1;   % change class label 0 -> -1
        y(y == 1) =  1;   % change class label 1 ->  1

        % extract predictors for both classes
        x1 = trainx(y == -1, 2:end)';  
        x2 = trainx(y ==  1, 2:end)';  

        % concatenate to form design matrix
        x = [x1 x2]'; 

        % number of samples in each class
        n1 = sum(y == -1); 
        n2 = sum(y == 1);  

        % get mean for both classes (m1 and m2) and all variables
        m1 = mean(x1, 2);
        m2 = mean(x2, 2);

        % calculate covariance matrix for both classes
        S1 = cov(x1');
        S2 = cov(x2');


        % define GP function for mean, covariance and likelihood    
        meanfunc = @meanConst; 

        % define corresponding hyper-parameters for these functions
        if strcmp(func2str(covfunc), 'covSEard')
            hyp.mean = 0;
            hyp.cov  = log(ones(1, nr_vars + 1));
        elseif strcmp(func2str(covfunc), 'covSEiso')
            hyp.mean = 0;
            hyp.cov  = log(ones(1, 2));  % length and variance parameter
        end



        % 3rd argument : limit of function evaluations
        hyp = minimize(hyp, @gp, -40, inffunc, meanfunc, covfunc, likfunc, x, y);

        [a b c d lp] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x, y, testx, 1);

        post_prob = exp(lp);

       
        fprintf('[LOOCV, sample %i] point t [%g %g] has p = %g\n', i, testx(1), testx(2), post_prob);
       
        class_probs = [class_probs, post_prob];


    end

    % sample 2-class Binomial given prediction probabilities
    if(strcmp(PREDICTION_METHOD, 'BINOM_SAMPLING'))
     
        predictions = binornd(1, class_probs);
    
    elseif (strcmp(PREDICTION_METHOD, 'DECISION_BOUNDARY'))
    
        % Use 0.5 decision boundary to decide the class 
        predictions = class_probs;
        predictions(predictions > 0.5) = 1;
        predictions(predictions ~= 1) = 0;
        
    end
        
    
    
    
    spec = Specificity(true_classes', predictions);
    sens = Sensitivity(true_classes', predictions);
    
    fprintf('Specificity: %g, Sensitivity: %g\n', spec, sens);
    
    % write to file if filename is not NULL
    if ~isempty(outfile_csv)
        dlmwrite(outfile_csv, predictions);
    end

end


% Sensitivity: 
%
% If a person has a disease, how often will the test be positive (true positive rate)? 
%
% Put another way, if the test is highly sensitive and the test result is negative you can be nearly certain that they don’t have disease. 
% A Sensitive test helps rule out disease (when the result is negative). Sensitivity rule out or "Snout"
%
% Sensitivity= true positives/(true positive + false negative)
%
function [res] = Sensitivity(true, test) 

    
    % number of true positives
    TP = sum(true & test);

    % number of false negatives
    FN = sum(true & ~test);

    res = TP / (TP + FN);
end


% Specificity: 
%
% If a person does not have the disease how often will the test be negative (true negative rate)?
%
% In other terms, if the test result for a highly specific test is positive you can be nearly certain that they actually have the disease.
% A very specific test rules in disease with a high degree of confidence Specificity rule in or "Spin".
%
% Specificity=true negatives/(true negative + false positives)
%
function [res] = Specificity(true, test) 

    % number of true negatives
    TN = sum(~true & ~test);

    % number of false negatives
    FP = sum(~true & test);

    res = TN / (TN + FP);
end


%
% Load data for different setups and transform for direct use in gp_core.
%
function [heart_data, vnames] = load_data(datafile, DATA_SETUP_TEST)

    %
    % Load data and transform
    %
    
    datatmp = importdata(datafile, ';', 1);
    all_data = datatmp.data;
    all_names = datatmp.colheaders;

    % Define column index for some variable names
    % Edit this for own names.
    id_type = find(ismember(all_names, 'Type')); 
    id_sbp = find(ismember(all_names, 'SBP'));
    id_ta = find(ismember(all_names, 'Ta'));
    id_treq = find(ismember(all_names, 'Treq'));
    id_edv = find(ismember(all_names, 'EDV'));
    id_esv = find(ismember(all_names, 'ES'));
    id_tatreq = find(ismember(all_names, 'Ta_Treq'));
    id_tasbp = find(ismember(all_names, 'Ta_SBP'));
    id_age   = find(ismember(all_names, 'Age'));
   
    heart_data = [];
    
    %
    % Extract columns for each data setup and put the variable names into a
    % string vector for later plotting.
    %
    if DATA_SETUP_TEST == 1
    
        heart_data = all_data(:,[id_type, id_ta, id_treq, id_sbp, id_edv, id_tatreq, id_tasbp]);
        vnames = {all_names{id_ta}, all_names{id_treq}, all_names{id_sbp}, all_names{id_edv}, all_names{id_tatreq}, all_names{id_tasbp} };
    
    elseif DATA_SETUP_TEST == 2
    
        heart_data = all_data(:,[id_type, id_ta, id_treq, id_sbp, id_edv]);
        vnames = {all_names{id_ta}, all_names{id_treq}, all_names{id_sbp}, all_names{id_edv}};
    
    elseif DATA_SETUP_TEST == 3
        
        heart_data = all_data(:,[id_type, id_ta, id_edv, id_tatreq, id_tasbp]);
        vnames = {all_names{id_ta}, all_names{id_edv}, all_names{id_tatreq}, all_names{id_tasbp} };
    
    elseif DATA_SETUP_TEST == 4  % Tried to add 'Age' with some difficulties..
        %
        % Tried to add 'Age' to the data set. It will not run if id_tasbp
        % is also in the set (!!). Why? I don't know yet. The Cholesky 
        % decomposition fails. Matrix is not invertable.
        %
        % THIS WILL FAIL ->
        %heart_data = all_data(:,[id_type, id_ta, id_treq, id_sbp, id_edv, id_tatreq, id_tasbp, id_age]);
        %vnames = {all_names{id_ta}, all_names{id_treq}, all_names{id_sbp}, all_names{id_edv}, all_names{id_tatreq}, all_names{id_tasbp}, all_names{id_age} };
        
        % But not this: (missing id_tasbp)
        heart_data = all_data(:,[id_type, id_ta, id_treq, id_sbp, id_edv, id_tatreq, id_age]);
        vnames = {all_names{id_ta}, all_names{id_treq}, all_names{id_sbp}, all_names{id_edv}, all_names{id_tatreq}, all_names{id_age} };
        
    end
    
    % Change mean to zero and std-dev to 1
    %  This is to make finding good hyper-parameters easier. 
    zscore_vars = bsxfun(@minus, heart_data(:,2:end), mean(heart_data(:,2:end)));
    zscore_vars = bsxfun(@rdivide, zscore_vars, std(zscore_vars));
    heart_data(:,2:end) = zscore_vars;



end


%
% This script creates a posterior probability plot of the Gaussian Process
% given the two most relevant input variables determined by Automatic
% Relevance Determination. 
%
% It also creates a plot for the different variables and their
% corresponding length scales using ARD.
%
%function gpard_lscales_probmesh(DATA_SETUP_TEST, datafile, OUT_DIR, PLOT_LSCALE_BARPLOT, likfunc)
function gpard_lscales_probmesh(heart_data, vnames, OUT_DIR, PLOT_LSCALE_BARPLOT, likfunc, DATA_SETUP_TEST)  
    
    % extract response and change labels
    y = heart_data(:,1);   % class response
    y(y == 0) = -1;   % change class label 0 -> -1
    y(y == 1) =  1;   % change class label 1 ->  1

    % extract predictors for both classes
    x1 = heart_data(y == -1, 2:end)';  
    x2 = heart_data(y ==  1, 2:end)';  

    % concatenate to form design matrix
    x = [x1 x2]'; 

    
    
    %
    % Run GP ARD and extract length scale parameters
    %
    
    
    meanfunc = @meanConst; 
    covfunc = @covSEard;   

    hyp.mean = 0;
    hyp.cov = log(ones(1, (size(x, 2) + 1) ));

    % minimize ARD length scale parameters
    hyp = minimize(hyp, @gp, -40, @infEP, meanfunc, covfunc, likfunc, x, y);

    % retrieve length-scales
    length_scales = hyp.cov(1:end-1);
    
    % sort length scales in ascending order. As smaller as more important 
    % is the corresponding parameter
    [~, sidx] = sort(length_scales);
        
    select_var1 = sidx(1);
    select_var2 = sidx(2);
    select_var3 = sidx(3);
    
    fprintf('The two most relevant variables are %s and %s\n', vnames{select_var1}, vnames{select_var2});
    
    
    %
    % Plot length scales relative to input range
    %

    relative_lscales = (length_scales ./ range(x));

    % invert it so it can be compared better to Lasso
    relative_lscales = 1 - relative_lscales; 

    if PLOT_LSCALE_BARPLOT
        % start plotting the bar plot
        figure(1);
        hold off;

        bar(relative_lscales,'DisplayName','normlscales');

        set(gca, 'FontSize', 16)
        title('GP ARD: inverse length scales')
        ylabel('importance (1 - relative length)');
        xlabel('variables');


        set(gca, 'FontSize', 10)
        set(gca,'XTickLabel', vnames);

        % save as pdf
        fileout = sprintf('%s/gpard_lengthscales_dataset%i.pdf', OUT_DIR, DATA_SETUP_TEST);
        set(gcf, 'PaperPosition', [0 0 12 12]); %Position plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [12 12]); %Set the paper to have width 5 and height 5.
        saveas(gcf, fileout, 'pdf')   %Save figure
    end
    
    % save numeric values in order to make consistent looking barplot with
    % the R script run_MI_classifiers.R
    head = vnames; 
    head=sprintf('%s,',head{:});
    head(end)='';
    
    % write the specificities to a .csv file
    outfile = sprintf('%s/gpard_lengthscales_dataset%i.csv', OUT_DIR, DATA_SETUP_TEST);
    dlmwrite(outfile, head, '');
    dlmwrite(outfile, relative_lscales, '-append')
    
    
    %
    % Plot the posterior probability mesh of the most important variable against
    % the second and third most important.
    %
      
    gpard_probmesh(x(:,[select_var1, select_var2]), y, vnames([select_var1, select_var2]), OUT_DIR, DATA_SETUP_TEST, likfunc);
    gpard_probmesh(x(:,[select_var1, select_var3]), y, vnames([select_var1, select_var3]), OUT_DIR, DATA_SETUP_TEST, likfunc);
    
    
    
end



%
% Plots a mesh grid of posterior probabilities from a Gaussian Process with 
% Automatic Relevance Determination (ARD) covariance function given 
% a 2-dimensional input x and two class output y.
%
%

function gpard_probmesh(x, y, labels, OUT_DIR, DATA_SETUP_TEST, likfunc)

    grid_points = 50;
    
    % boundaries for plot and probability mesh grid
    minx = min(x(:,1));
    maxx = max(x(:,1));
    miny = min(x(:,2));
    maxy = max(x(:,2));

    x1 = x(y == -1,:)';
    x2 = x(y ==  1,:)';
    
    % create probability mesh
    x_grid_interval = (maxx - minx) / grid_points;
    y_grid_interval = (maxy - miny) / grid_points;
    
    [t1 t2] = meshgrid(minx:x_grid_interval:maxx, miny:y_grid_interval:maxy);
    t = [t1(:) t2(:)]; 
    n = length(t);               % these are the test inputs

    meanfunc = @meanConst; 
    hyp.mean = 0;
    covfunc = @covSEard;   
    hyp.cov = log(ones(1, (size(x, 2) + 1) ));
    

    hyp = minimize(hyp, @gp, -40, @infEP, meanfunc, covfunc, likfunc, x, y);
    [a b c d lp] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, t, ones(n,1));

    % start plotting
    figure(2);
    hold off;
    
    set(gca, 'FontSize', 18)
    
    set(gcf, 'colormap', gray)
   
    % plot the probability mesh (contout)
    contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9, 'ShowText', 'on')
    hold on
    [~, h] = contour(t1, t2, reshape(exp(lp), size(t1)), [0.5 0.5], 'LineWidth', 2, 'ShowText', 'on');
    
    %colorbar
    
    grid
    axis([minx - 0.1 maxx + 0.1 miny - 0.05 maxy + 0.05])

    % plot points for each class 
    plot(x1(1,:), x1(2,:), '+', 'MarkerSize', 12, 'LineWidth',2, 'Color', [0.1 0.1 0.1]); 
    hold on
    plot(x2(1,:), x2(2,:), 'o', 'MarkerSize', 12, 'LineWidth',2, 'Color', [0.1 0.1 0.1]);

    title('GP ARD: posterior probabilities', 'FontSize', 16);

    % Change label format for feature names with subscript
    lx = labels(1);
    ly = labels(2);
    
    % .. for the xlabel
    lx = strrep(lx, '_', '/');
    lx = strrep(lx, 'Ta/SBP', 'Ta^{norm}');
    lx = strrep(lx, 'Ta/Treq', 'C^s'); %'Ta^{req}');
    lx = strrep(lx, 'Ta', 'T_a');
    lx = strrep(lx, 'Treq', 'T_{req}');
    
    % .. for the ylabel
    ly = strrep(ly, '_', '/');
    ly = strrep(ly, 'Ta/SBP', 'Ta^{norm}');
    ly = strrep(ly, 'Ta/Treq', 'C^s'); %'Ta^{req}');
    ly = strrep(ly, 'Ta', 'T_a');
    ly = strrep(ly, 'Treq', 'T_{req}');
    
    
    
    % set the labels
    xlabel(lx, 'FontSize', 14);
    ylabel(ly, 'FontSize', 14);
    
    % For the file name use yet another label, i.e. '_' instead of '/'
    format_label1 = labels{1};
    format_label2 = labels{2};
      
    % save as eps
    %fileout = sprintf('Results/GP_ARD/gpard_meshprob_%s_%s.eps', format_label1, format_label2);
    %print('-depsc', fileout);

    % save as pdf
    fileout = sprintf('%s/gpard_meshprob_%s_%s_dataset%i.pdf', OUT_DIR, format_label1, format_label2, DATA_SETUP_TEST);
    set(gcf, 'PaperPosition', [0 0 14 12]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [14 12]); %Set the paper to have width 5 and height 5.
    saveas(gcf, fileout, 'pdf')   %Save figure

    
end

