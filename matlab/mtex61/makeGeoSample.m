function sample = makeGeoSample(strike,dip,trend,varargin)
%MAKEGEOSAMPLE Keep field fabric and optional EBSD data in paired frames.
% Optional arguments are a plottingConvention, EBSD objects, and phase
% names. EBSD orientations are rotated while their map XY positions remain.

how2plot = plottingConvention.default;
ebsdInputs = {};
phaseNames = {};
for k = 1:numel(varargin)
    value = varargin{k};
    if isa(value,'plottingConvention')
        how2plot = value;
    elseif isa(value,'EBSD')
        ebsdInputs{end+1} = value; %#ok<AGROW>
    elseif ischar(value) || (isstring(value) && isscalar(value))
        phaseNames{end+1} = char(value); %#ok<AGROW>
    else
        error('makeGeoSample:Argument','Unsupported optional argument at position %d.',k);
    end
end

[fabOr,fabRot,strikeV,poleV,linV] = SDT2or(strike,dip,trend,how2plot);
[lineTrend,plunge] = V2TP(linV,how2plot);

sample.geo.strike = strike;
sample.geo.dip = dip;
sample.geo.trend = trend;
sample.geo.lineTrend = lineTrend;
sample.geo.plunge = plunge;
sample.geo.rake = angle(linV,strikeV,'noSymmetry')./degree;
sample.geo.pole = poleV;
sample.geo.lin = linV;
sample.geo.fabOr = fabOr;
sample.geo.fabRot = fabRot;
sample.geo.how2plot = how2plot;

if dot(how2plot.outOfScreen,zvector) > 0
    sample.spec.lin = xvector;
    sample.spec.pole = -yvector;
else
    sample.spec.lin = -xvector;
    sample.spec.pole = yvector;
end
sample.spec.how2plot = how2plot;

if ~isempty(phaseNames) && isempty(ebsdInputs)
    error('makeGeoSample:Phase','Phase names require an EBSD input.');
end
if ~isempty(ebsdInputs)
    sample.spec.ebsd = ebsdInputs;
    sample.geo.rotebsd = cell(size(ebsdInputs));
    for k = 1:numel(ebsdInputs)
        sample.geo.rotebsd{k} = rotEBSD2Geo(ebsdInputs{k},fabRot);
    end
    for k = 1:numel(phaseNames)
        phase = phaseNames{k};
        sample.geo.(phase) = rotTex2Geo(ebsdInputs{1}(phase).orientations,fabRot);
    end
end
end
