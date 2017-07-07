
%
%
% This script creates a posterior probability plot of the Gaussian Process
% given the two most relevant input variables determined by Automatic
% Relevance Determination. 
%
% It also creates a plot for the different variables and their
% corresponding length scales using ARD.
%

function gpard_lscales_prob() 

    clear all, close all

    DATA_SETUP_TEST = 1;
    write_figs = 1;

    % load the gpml package
    currpath = pwd();
    cd gpml-matlab;
    startup;
    cd(currpath);
    
    %
    % Load data and transform
    %
    
    datatmp = importdata('heart_ver2.csv', '\t', 1);
    all_x = datatmp.data;
    
    % get the names by using the old split method
    % ! Change the second argument if it is not a tabular separated table
    vnames = regexp(datatmp.textdata, '\t', 'split');
    vnames = vnames{1};    
        
    % the last element in then names could be empty, remove it
    vnames = vnames(1:(end-1));
    
    % add two new variables, ratios of the existing ones
    if DATA_SETUP_TEST == 1
        all_x = [all_x, (all_x(:,4) ./ all_x(:,3))];  % Ta/Tref
        all_x = [all_x, (all_x(:,4) ./ all_x(:,2))];  % Ta/SBP

        % append new names
        vnames{end+1} = 'Ta/Tref';
        vnames{end+1} = 'Ta/SBP';
    end
    
    %
    % Change mean to zero and std-dev to 1
    %  This is to make finding good hyper-parameters easier. 
    %
    zscore_vars = bsxfun(@minus, all_x(:,2:end), mean(all_x(:,2:end)));
    zscore_vars = bsxfun(@rdivide, zscore_vars, std(zscore_vars));
    all_x(:,2:end) = zscore_vars;
   
    % extract response and change labels
    y = all_x(:,1);   % class response
    y(y == 0) = -1;   % change class label 0 -> -1
    y(y == 1) =  1;   % change class label 1 ->  1

    % extract predictors for both classes
    x1 = all_x(y == -1, 2:end)';  
    x2 = all_x(y ==  1, 2:end)';  

    % concatenate to form design matrix
    x = [x1 x2]'; 

    
    
    %
    % Run GP ARD and extract length scale parameters
    %
    
    
    meanfunc = @meanConst; 
    covfunc = @covSEard;   
    likfunc = @likErf;

    hyp.mean = 0;
    hyp.cov = log(ones(1, (size(x, 2) + 1) ));

    % minimize ARD length scale parameters
    hyp = minimize(hyp, @gp, -40, @infEP, meanfunc, covfunc, likfunc, x, y);

    % retrieve length-scales
    length_scales = hyp.cov(1:end-1);
    
    % the two smallest length-scales matrch the two most relevant variables
    [~, sidx] = sort(length_scales);
        
    select_var1 = sidx(1);
    select_var2 = sidx(2);
    
    fprintf('The two most relevant variables are %s and %s\n', vnames{select_var1+1}, vnames{select_var2+1});
    
    
    %
    % Plot length scales relative to input range
    %
    
    relative_lscales = length_scales ./ range(x);
    
    figure(1);
    hold off;
    
    
    bar(relative_lscales,'DisplayName','normlscales');
    
    set(gca, 'FontSize', 16)
    title('GP ARD: length scales')
    ylabel('relative length');
    xlabel('variables');
    
       
    set(gca, 'FontSize', 10)
    set(gca,'XTickLabel', vnames(2:end));
    
    
    if write_figs
        % save as eps
        %fileout = sprintf('0_PLOTS/gpard_lengthscales.eps');
        %print('-depsc', fileout);

        
        % save as pdf
        fileout = sprintf('0_PLOTS/gpard_lengthscales_datatest%i.pdf', DATA_SETUP_TEST);
        
        set(gcf, 'PaperPosition', [0 0 12 12]); %Position plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [12 12]); %Set the paper to have width 5 and height 5.
        saveas(gcf, fileout, 'pdf')   %Save figure
        
    end

    
    
    
    %
    % Plot the posterior probability mesh with these two variables
    %
      
    %% TODO: VARIABLES SHOULD BE IN OLD SCALE AND MEAN. POST-PROCESS POINTS TO PLOT (TIMES OLD_SD, + OLD_MEAN) ?
    gpard_probmesh(x(:,[select_var1, select_var2]), y, vnames([select_var1+1, select_var2+1]), write_figs);
    
    
    
    
end