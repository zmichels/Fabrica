function [trd,plg] = V2TP(vec3d,how2plot)
%V2TP Convert MTEX vectors to clockwise trend and downward plunge.
% Angles are degrees. Axial vectors are brought into the screen's lower
% hemisphere before conversion. Horizontal directions retain their sign.

if nargin < 2 || isempty(how2plot)
    how2plot = plottingConvention.default;
end
if ~isa(how2plot,'plottingConvention')
    error('V2TP:Convention','Expected an MTEX plottingConvention.');
end

v = vec3d;
down = dot(v,how2plot.intoScreen);
len = sqrt(v.x.^2 + v.y.^2 + v.z.^2);
if any(len(:) == 0)
    error('V2TP:ZeroVector','Trend and plunge are undefined for zero vectors.');
end
v(down < -1e-10.*len) = -v(down < -1e-10.*len);
east = dot(v,how2plot.east);
north = dot(v,how2plot.north);
trd = mod(atan2d(east,north),360);
trd(hypot(east,north) < 1e-10.*len) = NaN;
plg = asind(max(-1,min(1,dot(v,how2plot.intoScreen)./len)));
trd = trd(:);
plg = plg(:);
end
