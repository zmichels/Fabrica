function [sample] = makeGeoSample(strike, dip, trend, varargin)

[fabOr, fabRot, strikeV, poleV, linV] = SDT2or(strike,dip,trend);

plunge = round(angle(linV,vector3d(linV.x,linV.y,0),'antipodal')./degree);

sample.geo.strike = strike;
sample.geo.dip = dip;
sample.geo.trend = trend;
sample.geo.plunge = plunge;
sample.geo.rake = round(angle(linV,strikeV,'antipodal')./degree);
sample.geo.pole = poleV;
sample.geo.lin = linV;
sample.geo.fabOr = fabOr;
sample.geo.fabRot = fabRot;

sample.spec.lin = xvector;
sample.spec.pole = -yvector;

numvars = size(varargin,2);
for n = 1:numvars
vclass{n} = class(varargin{n});
end

IEBSD = find(contains(vclass,'EBSD'));
ICHAR = find(contains(vclass,'char'));

if isempty(IEBSD)==0
    
    for j = 1:length(IEBSD)
        
        specEBSD{j,:} = varargin{IEBSD(j)};
        geoEBSD{j,:} = rotEBSD2Geo(specEBSD{j,:},fabOr);
        
    end
    
    sample.spec.ebsd = specEBSD;
    sample.geo.rotebsd = geoEBSD;
    
end
    


if isempty(ICHAR)==0
    
    for j = 1:length(ICHAR)
        
        poi{j} = varargin{ICHAR(j)};

        sample.geo.(poi{j}) = rotTex2Geo(sample.spec.ebsd{1}(poi{j}).orientations,fabOr);
       
        
    end
    

    
end