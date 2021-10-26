function plot_rec02(figNr, x, c, X, params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %%% Parse input arguments %%%
    switch nargin
        case 3
            X = [];
            params = {};

        case 4
            params = {};
        
        case 5

        otherwise
            error('Wrong number of inputs!')
    end
    %%%%%%%%%%%%%%%%%%%
    
    nrSteps = size(x,1) - 1;
    figure(figNr)
    colormap winter;
    scatter(x(1,:)*180/pi, x(2,:)*180/pi, 25, c, 'filled');
    if ~isempty(X)
        hold on
        X = X*(180/pi);
        X.plot('wire', true, 'linestyle', '-', 'linewidth', 2)
    end
    hold on;
    cmap = winter;
    p1 = scatter(NaN, NaN, 5, cmap(1,:), 'filled');
    p2 = scatter(NaN, NaN, 5, cmap(end,:), 'filled');
    xlabel('position [deg]')
    ylabel('velocity [deg/s]')
    legend([p1, p2], {'infeasible initial conditions', 'feasible initial conditions'}, 'Position',[0.75 0.91 0.08 0.06])
    set(gcf,'position',[100,100,params.width,1.5*params.height],'color','white')
end

