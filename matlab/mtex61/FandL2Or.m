function [fabOr, strikes, poles, lins] = FandL2Or(dipDirection,dip,rake,how2plot)
%FANDL2OR Fabric from dip direction, dip, and right-hand rake (degrees).
% Rake 0 follows strike; rake 90 follows the downward dip direction.
% This entry point preserves its historical specimen X=lineation,
% Z=foliation-pole reference axes.

if nargin < 4 || isempty(how2plot)
    how2plot = plottingConvention.default;
end
[~,~,strikes,poles,lins] = SDP2or(dipDirection-90,dip,rake,how2plot);
fabRot = rotation.map(zvector,poles,xvector,lins);
fabOr = orientation(fabRot,crystalSymmetry('mmm'));
end
