function [fabOr, fabRot, strikeV, poleV, linV] = geoFab(strike,dip,trend,linDir,poleDir,how2plot)
%GEOFAB Map explicit specimen lineation and foliation-pole axes to field axes.
% Strike is right-hand-rule; dip and trend are degrees. The lineation is
% the down-going intersection of foliation and the vertical trend plane.

if nargin < 6 || isempty(how2plot)
    how2plot = plottingConvention.default;
end
if any(dip(:) < 0 | dip(:) > 90)
    error('geoFab:Dip','Dip must be between 0 and 90 degrees.');
end
strikeV = TP2vec(strike,zeros(size(strike)),how2plot);
trendV = TP2vec(trend,zeros(size(trend)),how2plot);
poleV = TP2vec(strike-90,90-dip,how2plot);

% Keep the measured horizontal trend while solving its vertical component.
poleDown = dot(poleV,how2plot.intoScreen);
horizontalDot = dot(trendV,poleV);
verticalPlane = abs(poleDown) < 1e-10;
if any(verticalPlane(:) & abs(horizontalDot(:)) > 1e-10)
    error('geoFab:Trend','This trend is not in the vertical foliation plane.');
end
verticalPart = zeros(size(poleDown));
verticalPart(~verticalPlane) = -horizontalDot(~verticalPlane) ./ poleDown(~verticalPlane);
linV = normalize(trendV + verticalPart.*how2plot.intoScreen);

% An axial lineation is represented by its lower-hemisphere direction.
up = dot(linV,how2plot.intoScreen) < -1e-10;
linV(up) = -linV(up);

fabRot = rotation.map(linDir,linV,poleDir,poleV);
fabOr = orientation(fabRot,crystalSymmetry('mmm'));
end
