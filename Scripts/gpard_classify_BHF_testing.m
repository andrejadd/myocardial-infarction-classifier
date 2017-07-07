%
% This script is for testing the gpml GP classifier with ARD on the BHF data. 
% It is adapted from the demoClassification.m script provided by the gpml documentation:
%
% http://www.gaussianprocess.org/gpml/code/matlab/doc/ 
%
%


% load the GPstuff package
currpath = pwd();
cd gpml-matlab;
startup;
cd(currpath);
   

clear all, close all

write_fig = 1;

grid_points = 50;

%
% Select data version
%
%  1 : old data (27 healty cases, 6  MI patients)
%  2 : old data (27 healty cases, 11 MI patients)
%
% Note: for version 1, the covariance matrix for the MI patient data is
% going to be singular and predictions will fail in most cases. 
%
% Fix: If necessary (though version 1 is not used), make S2 non-singular by
% adding a small element to the diagonal.
%
data_version = 2;

%
% The type of classification to use: 
%   run_type = 1 directly calculates probabilities with given mean and
%                covariance. 
%   run_type = 2 is the standard ARD type (* I use this)
%   run_type = 3 uses FITC for large scale approximate inference. 
%                * this one makes use of inducing points 'u'.                       
% 
%
run_type = 3;


% read in data
if data_version == 1
    % load version 1 data
    datatmp = importdata('heart_ver1.csv');
    all_x = datatmp.data;
    
elseif data_version == 2
    % load version 2 data
    datatmp = importdata('heart_ver2.csv');
    all_x = datatmp.data;

end



% extract data for class 
y = all_x(:,1);   % class response
y(y == 0) = -1;   % change class label 0 -> -1
y(y == 1) =  1;   % change class label 1 ->  1

n1 = sum(y == -1); % number of class 1
n2 = sum(y == 1);  % number of class 2

% input variables
x = all_x(:,2:3);

% Change mean to zero and std-dev to 1
%  This is to make finding good hyper-parameters easier. 
zscore_vars = bsxfun(@minus,x, mean(x));
zscore_vars = bsxfun(@rdivide, zscore_vars, std(zscore_vars));
x = zscore_vars;

% extract predictors for both classes
x1 = x(y == -1, :)';  
x2 = x(y ==  1, :)';  

% get mean for both variables (columns in x1 and x2) and both classes (m1
% and m2)
m1 = mean(x1, 2);
m2 = mean(x2, 2);

% calculate covariance matrix for both classes
S1 = cov(x1');
S2 = cov(x2');

%
% DEBUG, leave-one-out
%
%x = x(2:end,:);
%y = y(2:end);
%n1 = n1 - 1;


% boundaries for plot and probability mesh grid
minx = min(x(:,1));
maxx = max(x(:,1));
miny = min(x(:,2));
maxy = max(x(:,2));




%
% Use type 1 classification to get posterior probability grid 
%
if run_type == 1
    
    
    figure(6)
    plot(x1(1,:), x1(2,:), 'b+', 'MarkerSize', 12); hold on
    plot(x2(1,:), x2(2,:), 'r+', 'MarkerSize', 12);


    % create probability mesh
    x_grid_interval = (maxx - minx) / grid_points;
    y_grid_interval = (maxy - miny) / grid_points;
    
    % create probability mesh
    [t1 t2] = meshgrid(minx:x_grid_interval:maxx, miny:y_grid_interval:maxy);
    t = [t1(:) t2(:)]; 
    n = length(t);               % these are the test inputs

    disp('tmm = bsxfun(@minus, t, m1'');')
    tmm = bsxfun(@minus, t, m1');

    disp('p1 = n1*exp(-sum(tmm*inv(S1).*tmm/2,2))/sqrt(det(S1));')
    p1 = n1*exp(-sum(tmm*inv(S1).*tmm/2,2))/sqrt(det(S1));

    disp('tmm = bsxfun(@minus, t, m2'');')
    tmm = bsxfun(@minus, t, m2');

    disp('p2 = n2*exp(-sum(tmm*inv(S2).*tmm/2,2))/sqrt(det(S2));')
    p2 = n2*exp(-sum(tmm*inv(S2).*tmm/2,2))/sqrt(det(S2));
    set(gca, 'FontSize', 24)

    disp('contour(t1, t2, reshape(p2./(p1+p2), size(t1)), [0.1:0.1:0.9])')
    contour(t1, t2, reshape(p2./(p1+p2), size(t1)), [0.1:0.1:0.9])

    [c h] = contour(t1, t2, reshape(p2./(p1+p2), size(t1)), [0.5 0.5]);
    set(h, 'LineWidth', 2)
    colorbar
    grid
    axis([minx-1 maxx+1 miny-0.05 maxy+0.05])
    fileout = sprintf('0_PLOTS/postprob_data%i_ctype%i', data_version, run_type );
    if write_fig, print('-depsc', fileout); end


    %
    % Make prediction for specific point
    %
    t = [127, 3.2];

    tmm = bsxfun(@minus, t, m1');
    p1 = n1*exp(-sum(tmm*inv(S1).*tmm/2,2))/sqrt(det(S1));

    tmm = bsxfun(@minus, t, m2');
    p2 = n2*exp(-sum(tmm*inv(S2).*tmm/2,2))/sqrt(det(S2));

    % class probability that t is in class 2 (1)
    post_prob = p2./(p1+p2);
    fprintf('point t [%g %g] has p=%g\n', t(1), t(2), post_prob);

end

%
% Use type 2 classification to get posterior probability grid
%

if run_type == 2
    MAKE_MESH_GRID = 1;
    if MAKE_MESH_GRID
        % create probability mesh
        [t1 t2] = meshgrid(minx:1:maxx,miny:0.1:maxy);
        t = [t1(:) t2(:)]; 
        n = length(t);               % these are the test inputs


        disp('meanfunc = @meanConst; hyp.mean = 0;')
        meanfunc = @meanConst; 
        hyp.mean = 100;
        disp('covfunc = @covSEard;   hyp.cov = log([1 1 1]);')
        covfunc = @covSEard;   hyp.cov = log([1 1 1]);
        disp('likfunc = @likErf;')
        likfunc = @likErf;
        disp(' ')

        disp('hyp = minimize(hyp, @gp, -40, @infEP, meanfunc, covfunc, likfunc, x, y);')
        hyp = minimize(hyp, @gp, -40, @infEP, meanfunc, covfunc, likfunc, x, y);

        disp('[a b c d lp] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, t, ones(n, 1));')
        [a b c d lp] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, t, ones(n,1));

        disp(' ')
        figure(7)
        set(gca, 'FontSize', 24)
        disp('plot(x1(1,:), x1(2,:), ''b+''); hold on')
        plot(x1(1,:), x1(2,:), 'b+', 'MarkerSize', 12); hold on
        disp('plot(x2(1,:), x2(2,:), ''r+'')')
        plot(x2(1,:), x2(2,:), 'r+', 'MarkerSize', 12)
        disp('contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9)')
        contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9)
        [c h] = contour(t1, t2, reshape(exp(lp), size(t1)), [0.5 0.5]);
        set(h, 'LineWidth', 2)
        colorbar
        grid
        axis([minx-1 maxx+1 miny-0.05 maxy+0.05])

        fileout = sprintf('0_PLOTS/postprob_data%i_ctype%i', data_version, run_type );
        if write_fig, print('-depsc', fileout); end
    
        if write_fig, print -depsc class_plot2.eps; end
    end
    
    %
    % Make prediction for specific point
    %
    t = [127, 3.2];

    meanfunc = @meanConst; 
    hyp.mean = 100;

    covfunc = @covSEard;   
    hyp.cov = log([1 1 1]);
    likfunc = @likErf;
    hyp = minimize(hyp, @gp, -40, @infEP, meanfunc, covfunc, likfunc, x, y);
    %[a b c d lp] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, t, ones(n,1));
    [a b c d lp] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, t, ones(1,1));

    post_prob = exp(lp);

    fprintf('point t [%g %g] has p=%g\n', t(1), t(2), post_prob);
end

if run_type == 3
    %
    % Use type 3 classification to get posterior probability grid
    %
    
    % generate 'u' points, this is a parameter to @covFITC
    disp('large scale classification using the FITC approximation')
    disp('[u1,u2] = meshgrid(linspace(-2,2,5)); u = [u1(:),u2(:)];')
    [u1,u2] = meshgrid(linspace(-2,2,5)); u = [u1(:),u2(:)]; clear u1; clear u2
    
    % create probability mesh
    x_grid_interval = (maxx - minx) / grid_points;
    y_grid_interval = (maxy - miny) / grid_points;
    
    % create probability mesh
    [t1 t2] = meshgrid(minx:x_grid_interval:maxx, miny:y_grid_interval:maxy);
    t = [t1(:) t2(:)]; 
    n = length(t);               % these are the test inputs

    
    %nu = size(u,1);
    % Define function and hyper-parameters for the GP
    meanfunc = @meanConst; 
    covfunc = @covSEard;   
    covfuncF = {@covFITC, {covfunc}, u};
    likfunc = @likErf;
    hyp.mean = 0;
    hyp.cov = log([1 1 1]);
    
    % one could also use @infFITC_Laplace
    inffunc = @infFITC_EP;                     
    
    
    hyp = minimize(hyp, @gp, -40, inffunc, meanfunc, covfuncF, likfunc, x, y);
    [a b c d lp] = gp(hyp, inffunc, meanfunc, covfuncF, likfunc, x, y, t, ones(n,1));
       
    
    % Plot input points and probability mesh
    figure(8)
    set(gca, 'FontSize', 20)
    plot(x1(1,:), x1(2,:), 'b+', 'MarkerSize', 12); hold on
    plot(x2(1,:), x2(2,:), 'r+', 'MarkerSize', 12)
   
    disp('contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9)')
    contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9)
    [c h] = contour(t1, t2, reshape(exp(lp), size(t1)), [0.5 0.5]);
    
    set(h, 'LineWidth', 2)
    disp('plot(u(:,1),u(:,2),''ko'', ''MarkerSize'', 12)')
    plot(u(:,1),u(:,2),'ko', 'MarkerSize', 12)
    colorbar
    grid

    axis([minx - 0.1 maxx + 0.1 miny - 0.05 maxy + 0.05])

    if write_fig, print -depsc class_plot3.eps; end
    
    
end