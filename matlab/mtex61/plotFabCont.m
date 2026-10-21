function f = plotFabCont(fabOr,how2plot,varargin)
%PLOTFABCONT Contour the foliation poles and lineations of fabric data.
% plotFabCont(fabOr,pc) retains FandL2Or's +X lineation/+Z pole defaults.
% plotFabCont(fabOr,pc,lineDir,poleDir) uses explicit specimen fabric axes.
% SDT2or/SDP2or use +X/-Y for +Z-out or -X/+Y for -Z-out conventions.

if nargin < 2 || isempty(how2plot)
    how2plot = plottingConvention.default;
end
[lineV,poleV] = resolveFabPlotVectors(fabOr,varargin{:});
f = figure;
plot(poleV,'antipodal','lower','contour','LineColor','k',how2plot)
hold on
plot(lineV,'antipodal','lower','contour','LineColor','r',how2plot)
end
