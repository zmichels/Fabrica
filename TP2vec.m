function [v] = TP2vec(trend,plunge)
%%

% query plotting preference/convention/setting:

% is xvector preference set to east or west?
xdir = getMTEXpref('xAxisDirection');
% is zvector preference set in or out of page?
zdir = getMTEXpref('zAxisDirection');


% adjust variables for different cases of settings

% x-axis preference cases:
switch xdir
    
    case 'east'
        % no need to do anything
        
    case 'west'
        plunge = -plunge;
        
    case 'north'
        trend = trend + 90;
        
    case 'south'
        trend = trend - 90;
end
    
 
% z-axis preference cases:
switch zdir
    
    case 'outOfPlane'
        % no need to do anything
        
    case 'intoPlane'
        trend = -trend;
        plunge = -plunge;
    
end


% compute a vector3d object
v = vector3d('polar',(90+plunge)*degree,(90-trend)*degree);

v(v.z>0) = -v(v.z>0);
