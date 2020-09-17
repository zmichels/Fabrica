function [v] = TP2vec(trend,plunge)

% is xvector preference set to east or west?
xdir = getMTEXpref('xAxisDirection');
if strcmp(xdir,'east') == 0
    plunge = -plunge;
end

zdir = getMTEXpref('zAxisDirection');
if strcmp(zdir,'outOfPlane') == 0
    trend = -trend;
    plunge = -plunge;
end

v = vector3d('polar',(90+plunge)*degree,(90-trend)*degree);

if v.z>0
    v = -v;
end
