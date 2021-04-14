%% an example

% First plot the plotting directions you use to be certain we map to the
% correct specimen directions:
figure,
plot([vector3d.X,vector3d.Y,vector3d.Z],[1, 2, 3])
mtexColorMap(cool);


%% Setup and use new geoFab() function

% this function is effectively identical to the SDT2or() function in the
% Fabrica toolbox, except that it does not presume a specimen orientation,
% and instead, you can indicate the directions of fabric elements
% (lineation and foliation pole) at the input as vector3d objects.

% Let's say you generate the plot above and you see that the X-direction
% shown matches the lineation direction in your sample... you might then
% define:

linDir = vector3d.X;

% And let's say for this fictional sample that the *downward* pole to
% folaition points "down" on the page, opposite the Y-direction... you
% could then define:

poleDir = -vector3d.Y;


% define strike, dip, trend
strike = 110;
dip = 70;
trend = 15;

% run the new function
[fabOr, fabRot, strikeV, poleV, linV] = geoFab(strike,dip,trend,linDir,poleDir);


%% In a different case...

% let's say the lineation points opposite the X-direction in the plot you
% made at the start... then you could define "linDir" as the negative of
% the X-direction
linDir = -vector3d.X;

% and maybe the foliation pole still points down in this sample, too...
poleDir = -vector3d.Y;

% run the function
[fabOr, fabRot, strikeV, poleV, linV] = geoFab(strike,dip,trend,linDir,poleDir);

% Additionally, this generic form will allow to define directions that
% are *not* aligned with the principal axes/directions X, Y, or Z.