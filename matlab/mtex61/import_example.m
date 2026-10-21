%% Read the bundled synthetic example, independent of the current folder.
here = fileparts(mfilename('fullpath'));
fabricPairs = readtable(fullfile(here,'example.csv'));

%% Parse the data into separate variables;

% Dip direction is presumed to be an azimuth angle in degree
dd = fabricPairs.dipDir;

% Dip is the foliation-plane inclination downward from horizontal in degrees
d = fabricPairs.dip;

% rake is presumed to be a clockwise angle in degrees down the foliation
% plane from horizontal
r = fabricPairs.rhrake;
