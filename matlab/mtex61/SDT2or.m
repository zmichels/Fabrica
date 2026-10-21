function [fabOr, fabRot, strikeV, poleV, linV] = SDT2or(strike,dip,trend,how2plot)
%SDT2OR Fabric from right-hand-rule strike, dip, and lineation trend.
% Angles are degrees. The lineation plunge is solved by intersecting its
% vertical trend plane with the foliation plane. By default the specimen
% lineation is +X and foliation pole is -Y when +Z is out of screen;
% both signs reverse when +Z points into screen.

if nargin < 4 || isempty(how2plot)
    how2plot = plottingConvention.default;
end
if dot(how2plot.outOfScreen,zvector) > 1-1e-10
    linDir = xvector;
    poleDir = -yvector;
elseif dot(how2plot.outOfScreen,-zvector) > 1-1e-10
    linDir = -xvector;
    poleDir = yvector;
else
    error('SDT2or:Convention','Specimen Z must be normal to the screen.');
end
[fabOr, fabRot, strikeV, poleV, linV] = ...
    geoFab(strike,dip,trend,linDir,poleDir,how2plot);
end
