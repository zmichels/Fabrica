function rodf = rotODFGeo(odf,fabric)
%ROTODFGEO Rotate an MTEX ODF using a fabric rotation or sample structure.
if isstruct(fabric)
    if ~isfield(fabric,'geo') || ~isfield(fabric.geo,'fabRot')
        error('rotODFGeo:Sample','Sample must contain geo.fabRot.');
    end
    fabric = fabric.geo.fabRot;
end
if ~isa(odf,'SO3Fun') || ~isa(fabric,'rotation')
    error('rotODFGeo:Type','Expected an ODF and a fabric rotation.');
end
rodf = rotate(odf,rotation(fabric));
end
