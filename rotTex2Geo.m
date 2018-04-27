function [rTex] = rotTex2Geo(o,fabOr)


%% Rotate orientation data from specimen to geographic reference frame
warning off;
% Rotate the orientations
rTex = fabOr.*o;
warning on