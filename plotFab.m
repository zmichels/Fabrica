function [f] = plotFab(fabOr)

% This function uses the MTEX toolbox plotting parameters to display an
% orientation with a great circle and a point. The axes of the orientation
% that are plotted are selected based on other functions in the FabOr
% library. This function will plot one or more fabric orientations in this
% manner.

% Dependencies: MTEX toolbox

%% make a plot
f = figure;
plot(fabOr*zvector,'antipodal','lower','plane')
hold on
plot(fabOr*xvector,'antipodal','lower','Marker','d','MarkerZise',12,'MarkerFaceColor','k')