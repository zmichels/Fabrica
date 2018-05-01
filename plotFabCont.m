function [f] = plotFabCont(fabOr)

% This function uses the MTEX toolbox plotting parameters to display a
% contour plot of foliation and lineation components of a tectonite fabric

% Dependencies: MTEX toolbox

%% make a plot
f = figure;
plot(fabOr*zvector,'antipodal','lower','contour','LineColor','k')
hold on
plot(fabOr*xvector,'antipodal','lower','contour','LineColor','r')