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

% is xvector preference set to east or west?
xdir = getMTEXpref('xAxisDirection');
% is zvector preference set in or out of page?
zdir = getMTEXpref('zAxisDirection');


switch zdir
    
    case 'outOfPlane'
        poleDir = -yvector;
        linDir = xvector;
        
    case 'intoPlane'
        poleDir = yvector;
        linDir = -xvector;
    
end

% Use fabric vectors to define principal axes of a fabric orientation
fabRot = rotation('map',linDir, linV, poleDir, poleV);
fabOr = orientation(fabRot,crystalSymmetry,'mmm');
    
