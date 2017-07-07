%
% Plots a mesh grid of posterior probabilities from a Gaussian Process with 
% Automatic Relevance Determination (ARD) covariance function given 
% a 2-dimensional input x and two class output y.
%
%

function gpard_probmesh(x, y, labels, write_fig)

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
    likfunc = @likErf;

    hyp = minimize(hyp, @gp, -40, @infEP, meanfunc, covfunc, likfunc, x, y);
    [a b c d lp] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, t, ones(n,1));

    figure(2);
    hold off;
    
    set(gca, 'FontSize', 18)
    
    plot(x1(1,:), x1(2,:), 'b+', 'MarkerSize', 12); hold on
    plot(x2(1,:), x2(2,:), 'r+', 'MarkerSize', 12)
    
    contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9)
    [c h] = contour(t1, t2, reshape(exp(lp), size(t1)), [0.5 0.5]);
    set(h, 'LineWidth', 2)
    colorbar
    grid
    axis([minx - 0.1 maxx + 0.1 miny - 0.05 maxy + 0.05])

    xlabel(labels(1));
    ylabel(labels(2));
    
    title('GP ARD: posterior probabilities');
    
    if write_fig
        format_label1 = strrep(labels{1}, '/', '_');
        format_label2 = strrep(labels{2}, '/', '_');

        % save as eps
        %fileout = sprintf('0_PLOTS/gpard_meshprob_%s_%s.eps', format_label1, format_label2);
        %print('-depsc', fileout);
        
        % save as pdf
        fileout = sprintf('0_PLOTS/gpard_meshprob_%s_%s.pdf', format_label1, format_label2);
        set(gcf, 'PaperPosition', [0 0 14 14]); %Position plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [14 14]); %Set the paper to have width 5 and height 5.
        saveas(gcf, fileout, 'pdf')   %Save figure
    end

    
end

