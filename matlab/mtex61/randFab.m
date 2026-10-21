function [s,d,dd,t,p,r] = randFab(num)
%RANDFAB Random right-hand-rule foliation/lineation field pairs (degrees).
% The rake r runs from right-hand strike toward down dip, 0 to 180 degrees.
% Trend t is clockwise from north and plunge p is downward. A vertical
% lineation has no defined trend (NaN). This numerical generator does not
% depend on a plottingConvention; TP2vec(t,p,pc) places its output in MTEX.

validateattributes(num,{'numeric'},{'scalar','integer','nonnegative'}, ...
    mfilename,'num');
s = 360*rand(num,1);
d = 90*rand(num,1);
r = 180*rand(num,1);
dd = mod(s+90,360);

% Components in geographic East/North/Down. atan2d retains the correct
% quadrant when the rake passes 90 degrees (opposite-strike lineations).
east = cosd(r).*sind(s) + sind(r).*cosd(d).*cosd(s);
north = cosd(r).*cosd(s) - sind(r).*cosd(d).*sind(s);
down = sind(r).*sind(d);
t = mod(atan2d(east,north),360);
p = atan2d(down,hypot(east,north));
t(hypot(east,north) < 1e-10) = NaN;
end
