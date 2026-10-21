function [lineV,poleV] = resolveFabPlotVectors(fabOr,varargin)
%RESOLVEFABPLOTVECTORS Map specimen fabric axes for the plot helpers.
% [lineV,poleV] = resolveFabPlotVectors(fabOr) retains the historical
% FandL2Or specimen frame: +X lineation and +Z foliation pole.
% [lineV,poleV] = resolveFabPlotVectors(fabOr,lineDir,poleDir) uses explicit
% specimen axes, e.g. +X/-Y for SDT2or or SDP2or when +Z is out of screen.

if nargin == 1
    lineDir = xvector;
    poleDir = zvector;
elseif nargin == 3
    lineDir = varargin{1};
    poleDir = varargin{2};
else
    error('resolveFabPlotVectors:Axes', ...
        'Supply both specimen lineation and pole axes, or neither.');
end
if ~isa(fabOr,'orientation') || ~isa(lineDir,'vector3d') || ...
        ~isa(poleDir,'vector3d')
    error('resolveFabPlotVectors:Type', ...
        'Expected an MTEX orientation and vector3d specimen axes.');
end
lineV = fabOr * lineDir;
poleV = fabOr * poleDir;
end
