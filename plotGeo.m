function [] = plotGeo(strike,dip,trend)

[~, ~, ~, poleV, linV] = SDT2or(strike,dip,trend);


f = figure;
plot(poleV,'antipodal','lower','plane')
hold on
plot(linV,'antipodal','lower','Marker','d','MarkerZise',12,'MarkerFaceColor','k')
hold on
% plot(poleV,'antipodal','lower','Marker','o','MarkerZise',12,'MarkerEdgeColor','k','MarkerFaceColor','none')