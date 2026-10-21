function [] = plotGeo(strike,dip,trend,how2plot)

if nargin < 4 || isempty(how2plot)
    how2plot = plottingConvention.default;
end
[~, ~, ~, poleV, linV] = SDT2or(strike,dip,trend,how2plot);


f = figure;
plot(poleV,'antipodal','lower','plane',how2plot)
hold on
plot(linV,'antipodal','lower','Marker','d','MarkerSize',12,'MarkerFaceColor','k',how2plot)
hold on
