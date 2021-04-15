function [trd,plg]=V2TP(vec3d)

%% A function to convert vector3d data from MTEX into a trend and plunge

% Get size and dimensions of inpout data
sz = size(vec3d);

% initialize output variable
trd = NaN(sz);

% polar coordinates
[theta, rho, ~] = polar(vec3d);
theta = (theta./degree);
rho = (rho./degree);
plg = 90 - (180 - theta);
trd = (360 - rho) + 90;

% Set trend between 0 - 360 degrees
trd(trd>360) = trd(trd>360)-360;
trd(trd<0) = trd(trd<0)+360;

end