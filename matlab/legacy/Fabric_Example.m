%% USE CASES:
% You have a sample collected from a field location and you measured the
% orientations of fabric defining elements such as foliation, and
% lineation. You have also collected EBSD data from that sample. You would
% like to rotate some results (i.e., crystal orientation data, ODF,
% grainboundary directions, etc.) into a geographic reference fram using
% the geographic orientation of the samlpe's fabric.

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

% IMPORTANT NOTE: currently, this set of scripts does not require plunge
% because it assumes that the lineation lies within the foliation plane
% defined by strike and dip. Therefore, the plunge and rake are solved by
% projection of the lineation trend onto the foliation plane. *Eventually*s
% this script will be able to handle any valid comination of inputs to
% define a 3D fabric.

%% Define an orientation and set of MTEX vectors to represent the fabric
% Use function SDT2or
[fabOr, fabRot, strikeV, poleV, linV] = SDT2or(strike,dip,trend);


%% Alternatively, create a 'sample' data structure store fabric info
% One advantage of the 'sample' data structure is that it allows for any
% associated EBSD data to be stored in a 'specimen' container within the
% structure. If no EBSD data are specified during the following function
% call, then no 'specimen' container will be created with the strucutre. If
% however, an EBSD dataset is included, the structure will include these
% data in the 'spec'  (specimen) container, *and* it will include a geographically
% rotated version of the data in the 'geo' (geographic) container. Please
% note than when EBSD data are rotated, only the crystal orientations are fully
% rotated (i.e., the EBSD map will not be properly rotated into the desired
% perspective). Additionally, if a string is included for a phase (or
% phases) of interest from the ebsd dataset, a geographically rotated
% version of those orientations will be included in the 'geo' container.

[sample] = makeGeoSample(strike, dip, trend, ebsd, 'fo');


%% Now what?!
% View the contents of variable 'sample'



%% Make a plot of the *original* (not rotated) forsterite data
o = sample.spec.ebsd{1}('fo').orientations;
oCS = o.CS;
ho = [Miller(1,0,0,'hkl',oCS),Miller(0,1,0,'hkl',oCS),Miller(0,0,1,'hkl',oCS)];

figure,
plotPDF(o,h,'antipodal','lower','contourf','colorrange','equal')
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
plotPDF(sample.geo.fo,h,'antipodal','lower','contourf','colorrange','equal')
hold on

% show the orientation of specimen directions assumed for fabric rotation
% lineation is assumed parallel to specimen x-direction vector
annotate(sample.geo.lin,'Marker','d')

% foliation plane is assumed normal to specimen y-direction vector
annotate(sample.geo.pole,'plane')

% pole points "downward"
annotate(sample.geo.pole,'Marker','s')

% IMPORTANT NOTE: In the future, these scripts will allow you to defined
% any vector for lineation and pole to foliation in specimen reference
% frame, rather than assume/require that lineation points towards the
% xvector, with the foliation pole pointing towards the -yvector.