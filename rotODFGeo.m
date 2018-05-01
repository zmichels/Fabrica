function [rodf] = rotODFGeo(odf,varargin)

numvars = size(varargin,2);
for n = 1:numvars
vclass{n} = class(varargin{n});
end

ISTRUCT = find(contains(vclass,'struct'));
IROT = find(contains(vclass,'rotation'));

if sum(ISTRUCT)==1 && sum(IROT)==0
    
    fabRot = varargin{ISTRUCT}.geo.fabOr;
    rodf = rotate(odf,fabRot);
    
elseif sum(IROT)==1 && sum(ISTRUCT)==0
    
    fabRot = varargin{IROT};
    rodf = rotate(odf,fabRot);
    
else
    
    fprintf('\n\n!!!  ERROR !!!\nThe second input variable is *not* an appropriate type\nPlease check that it is either a ''sample'' structure or a ''rotation'' class object\n')
    rodf = [];
    
end

