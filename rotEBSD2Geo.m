function [rEBSD] = rotEBSD2Geo(ebsd,fabOr)


%% Rotate orientation data from specimen to geographic reference frame
warning off;
% Rotate the orientations
rEBSD = fabOr.*ebsd;
warning on