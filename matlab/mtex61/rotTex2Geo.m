function rTex = rotTex2Geo(o,fabRot)
%ROTTEX2GEO Apply a specimen-to-geographic rotation to orientations.
if ~isa(o,'orientation') || ~isa(fabRot,'rotation')
    error('rotTex2Geo:Type','Expected orientations and an MTEX rotation.');
end
rTex = orientation(rotation(fabRot).*rotation(o),o.CS,o.SS);
end
