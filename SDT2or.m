function [fabOr, fabRot, strikeV, poleV, linV] = SDT2or(strike,dip,trend)

% This function takes three measurements of strike, trend, and plunge
% in degrees and uses them to compute a rotational representation of a rock
% fabric (i.e., a single quantity that combines foliation and lineation).

% Dependencies: MTEX toolbox, TP2vec

%% convert to vectors and compute poles to foliation
% fabric dataset of strike, dip, trend, plunge to MTEX vectors
strikeV = TP2vec(strike,zeros(size(strike)));
trendV = TP2vec(trend,zeros(size(trend)));
poleV = TP2vec(strike-90,90-dip);
linV = cross(cross(trendV,zvector),poleV);

poleV(poleV.z>0) = -poleV(poleV.z>0);
linV(linV.z>0) = -linV(linV.z>0);


%% Convert to orientation/rotation
% Use fabric vectors to define principal axes of a fabric orientation
fabRot = rotation('map',xvector, linV, zvector, poleV);
fabOr = orientation(fabRot,crystalSymmetry,'mmm');
    
