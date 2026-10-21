function f = plotFab(fabOr,how2plot,varargin)
%PLOTFAB Plot the foliation plane and lineation of a fabric orientation.
% plotFab(fabOr,pc) retains FandL2Or's +X lineation/+Z pole defaults.
% plotFab(fabOr,pc,lineDir,poleDir) uses the actual specimen fabric axes.
% For SDT2or/SDP2or, pass +X/-Y when pc.outOfScreen is +Z, or -X/+Y
% when it is -Z. Supply the same plottingConvention used to build fabOr.

if nargin < 2 || isempty(how2plot)
    how2plot = plottingConvention.default;
end
[lineV,poleV] = resolveFabPlotVectors(fabOr,varargin{:});
f = figure;
plot(poleV,'antipodal','lower','plane',how2plot)
hold on
plot(lineV,'antipodal','lower','Marker','d','MarkerSize',12, ...
    'MarkerFaceColor','k',how2plot)
end
