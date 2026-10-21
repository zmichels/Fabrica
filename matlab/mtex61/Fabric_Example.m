%% USE CASES:
% This example uses MTEX's forsterite EBSD dataset as a specimen. A measured
% foliation and lineation define a specimen-to-geographic rotation. The
% script keeps the EBSD map XY coordinates fixed while rotating orientations.

%% Load some ebsd data
% For this example, we will use the "Forsterite" dataset included with MTEX
mtexdata forsterite;
cs = ebsd('f').CS;
h = [Miller(1,0,0,'direction',cs),Miller(0,1,0,'direction',cs),Miller(0,0,1,'direction',cs)];
%% Enter some information about the sample's fabric from field collection
% foliation information
strike = 15
dip = 56

% lineation information
trend = 42

% SDT2or solves the plunge from the intersection of the foliation plane
% and the vertical plane containing the measured lineation trend.
pc = plottingConvention(zvector,xvector); % east +X, north +Y, down -Z

%% Define an orientation and set of MTEX vectors to represent the fabric
% Use function SDT2or
[fabOr, fabRot, strikeV, poleV, linV] = SDT2or(strike,dip,trend,pc);


%% Store source and geographic products together
% sample.spec records the specimen axes and source EBSD map. sample.geo
% retains the field directions, rotation, rotated EBSD orientations, and
% selected phase orientations. EBSD map XY positions remain unchanged.

[sample] = makeGeoSample(strike, dip, trend, ebsd, 'fo', pc);


%% Now what?!
% View the contents of variable 'sample'



%% Make a plot of the *original* (not rotated) forsterite data
o = sample.spec.ebsd{1}('fo').orientations;
figure,
plotPDF(o,h,'antipodal','lower','smooth','halfwidth',10*degree,'colorRange','equal',pc)
hold on

% show the orientation of specimen directions assumed for fabric rotation
% lineation is assumed parallel to specimen x-direction vector
annotate(sample.spec.lin,'Marker','d')

% foliation plane is assumed normal to specimen y-direction vector
annotate(sample.spec.pole,'plane')

% pole points "downward"
annotate(sample.spec.pole,'Marker','s')


%% Make a plot of the geographically rotated EBSD data showing fabric
figure,
plotPDF(sample.geo.fo,h,'antipodal','lower','smooth','halfwidth',10*degree,'colorRange','equal',pc)
hold on

% show the orientation of specimen directions assumed for fabric rotation
% lineation is assumed parallel to specimen x-direction vector
annotate(sample.geo.lin,'Marker','d')

% foliation plane is assumed normal to specimen y-direction vector
annotate(sample.geo.pole,'plane')

% pole points "downward"
annotate(sample.geo.pole,'Marker','s')

% For a different specimen lineation/pole pair, call geoFab with explicit
% specimen vectors. The supplied pair must be perpendicular.
