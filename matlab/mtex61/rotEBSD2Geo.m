function rEBSD = rotEBSD2Geo(ebsd,fabRot)
%ROTEBSD2GEO Rotate EBSD orientations but keep map XY positions unchanged.
if ~isa(ebsd,'EBSD') || ~isa(fabRot,'rotation')
    error('rotEBSD2Geo:Type','Expected EBSD data and an MTEX rotation.');
end
rEBSD = rotate(ebsd,rotation(fabRot),'keepXY');
end
