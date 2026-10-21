function v = TP2vec(trend,plunge,how2plot)
%TP2VEC Geographic trend and downward plunge as an MTEX vector3d.
% Angles are degrees; trend is clockwise from screen north. A positive
% plunge points into the screen's lower hemisphere. Supply how2plot (for
% example, ebsd.how2plot) to avoid depending on the session default.

if nargin < 3 || isempty(how2plot)
    how2plot = plottingConvention.default;
end
if ~isa(how2plot,'plottingConvention')
    error('TP2vec:Convention','Expected an MTEX plottingConvention.');
end
if any(abs(plunge(:)) > 90)
    error('TP2vec:Plunge','Plunge must be between -90 and 90 degrees.');
end

v = (cosd(plunge).*sind(trend)).*how2plot.east + ...
    (cosd(plunge).*cosd(trend)).*how2plot.north + ...
    sind(plunge).*how2plot.intoScreen;
end
