%% geoFab: explicit specimen fabric axes to geographic field directions
% Run after startup_mtex and addpath on the Fabrica folder. All angles are
% degrees. The right-hand-rule strike, dip and lineation trend describe a
% field fabric; linDir/poleDir describe the matching directions in the
% specimen frame. Unlike SDT2or, geoFab makes no fixed-axis assumption.

pc = plottingConvention(zvector,xvector); % east +X, north +Y, down -Z
strike = 110;
dip = 70;
trend = 15;

%% Example 1: lineation +X, downward foliation pole -Y in the specimen
linDir = xvector;
poleDir = -yvector;
[fabOr,fabRot,strikeV,poleV,linV] = ...
    geoFab(strike,dip,trend,linDir,poleDir,pc); %#ok<ASGLU>
mappedLine = fabRot*linDir;
mappedPole = fabRot*poleDir;
assert(norm([mappedLine.x-linV.x,mappedLine.y-linV.y,mappedLine.z-linV.z]) < 1e-10);
assert(norm([mappedPole.x-poleV.x,mappedPole.y-poleV.y,mappedPole.z-poleV.z]) < 1e-10);
assert(abs(dot(linV,poleV)) < 1e-10);

%% Example 2: same field fabric from a differently aligned specimen
linDir = -xvector;
poleDir = -yvector;
[otherOr,otherRot,~,otherPole,otherLine] = ...
    geoFab(strike,dip,trend,linDir,poleDir,pc); %#ok<ASGLU>
assert(norm([otherLine.x-linV.x,otherLine.y-linV.y,otherLine.z-linV.z]) < 1e-10);
assert(norm([otherPole.x-poleV.x,otherPole.y-poleV.y,otherPole.z-poleV.z]) < 1e-10);

% The two rotations differ in specimen coordinates, yet both mapped axes
% represent the same physical geographic lineation and foliation plane.
