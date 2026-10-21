function suffix = fabricaMtexArtifactSuffix
%FABRICAMTEXARTIFACTSUFFIX Keep published names while isolating later versions.
version = getMTEXpref('version');
switch version
    case 'mtex-6.1.1'
        suffix = 'mtex61';
    case 'mtex-7.2.0'
        suffix = 'mtex72';
    otherwise
        assert(startsWith(version,'mtex-'),'Unrecognized MTEX version: %s',version);
        suffix = ['mtex',regexprep(version(6:end),'[^A-Za-z0-9]','')];
end
end
