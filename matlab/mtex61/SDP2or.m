function [fabOr, fabRot, strikeV, poleV, linV] = SDP2or(strike,dip,pitch,how2plot)
%SDP2OR Fabric from right-hand-rule strike, dip, and rake from strike.
% Angles are degrees. Positive pitch rotates from strike about the down
% foliation pole. The source specimen axes match SDT2or.

if nargin < 4 || isempty(how2plot)
    how2plot = plottingConvention.default;
end
if any(dip(:) < 0 | dip(:) > 90)
    error('SDP2or:Dip','Dip must be between 0 and 90 degrees.');
end
strikeV = TP2vec(strike,zeros(size(strike)),how2plot);
poleV = TP2vec(strike-90,90-dip,how2plot);
linV = rotate(strikeV,rotation.byAxisAngle(poleV,pitch*degree));
up = dot(linV,how2plot.intoScreen) < -1e-10;
linV(up) = -linV(up);

if dot(how2plot.outOfScreen,zvector) > 1-1e-10
    linDir = xvector;
    poleDir = -yvector;
elseif dot(how2plot.outOfScreen,-zvector) > 1-1e-10
    linDir = -xvector;
    poleDir = yvector;
else
    error('SDP2or:Convention','Specimen Z must be normal to the screen.');
end
fabRot = rotation.map(linDir,linV,poleDir,poleV);
fabOr = orientation(fabRot,crystalSymmetry('mmm'));
end
