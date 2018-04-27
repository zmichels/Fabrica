function [fabOr, fabRot, strikeV, poleV, linV] = SDT2or(strike,trend,plunge)

% This function takes three measurements of strike, trend, and plunge
% in degrees and uses them to compute a rotational representation of a rock
% fabric (i.e., a single quantity that combines foliation and lineation).

% Dependencies: MTEX toolbox, TP2vec

%% convert to vectors and compute poles to foliation
% fabric dataset of strike, dip, trend, plunge to MTEX vectors
linV = TP2vec(trend,plunge);
strikeV = TP2vec(strike,zeros(size(strike)));
poleV = cross(linV,strikeV);
poleV(poleV.z>0) = -poleV(poleV.z>0);


%% Convert to orientation/rotation
% Use fabric vectors to define principal axes of a fabric orientation
fabRot = rotation('map',xvector, linV, -zvector, poleV);
fabOr = orientation(fabRot,crystalSymmetry,'mmm');
    
