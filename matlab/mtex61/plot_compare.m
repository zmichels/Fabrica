function [fSpec,fGeo] = plot_compare(o,rTex,h,fabOr,how2plot,varargin)
%PLOT_COMPARE Compare specimen and rotated pole figures with fabric axes.
% plot_compare(o,rTex,h,fabOr,pc) uses the historical +X line/+Z pole
% specimen frame of FandL2Or. To plot SDT2or/SDP2or results, supply both
% specimen axes: plot_compare(o,rTex,h,fabOr,pc,lineDir,poleDir).
% Use the same plottingConvention used to construct fabOr and rTex.

if nargin < 5 || isempty(how2plot)
    how2plot = plottingConvention.default;
end
if isempty(varargin)
    lineDir = xvector;
    poleDir = zvector;
elseif numel(varargin) == 2
    lineDir = varargin{1};
    poleDir = varargin{2};
else
    error('plot_compare:Axes', ...
        'Supply both specimen lineation and pole axes, or neither.');
end
[lineV,poleV] = resolveFabPlotVectors(fabOr,lineDir,poleDir);

fSpec = figure;
plotPDF(o,h,'antipodal','lower','smooth','colorRange','equal',how2plot)
hold on
annotate(lineDir,'antipodal','lower','Marker','d','MarkerSize',12, ...
    'MarkerFaceColor','k','add2all')
annotate(poleDir,'antipodal','lower','plane','add2all')
hold off

fGeo = figure;
plotPDF(rTex,h,'antipodal','lower','smooth','colorRange','equal',how2plot)
hold on
annotate(lineV,'antipodal','lower','Marker','d','MarkerSize',12, ...
    'MarkerFaceColor','k','add2all')
annotate(poleV,'antipodal','lower','plane','add2all')
hold off
end
