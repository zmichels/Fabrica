function [f] = plotFabCont(fab)

% This function uses the MTEX toolbox plotting parameters to display a
% contour plot of foliation and lineation components of a tectonite fabric

% Dependencies: MTEX toolbox

%% make a plot
f = figure;
plot(fab*zvector,'antipodal','lower','contour','LineColor','k')
hold on
plot(fab*xvector,'antipodal','lower','contour','LineColor','r')