function [rTex] = rotTex2Geo(o,fabOr)

%% Rotate orientation data from specimen to geographic reference frame
warning off;

% Check that the input is an orientation
check1 = strcmp(class(o),'orientation');
check2 = strcmp(class(fabOr),'orientation');

if check1 == 1 && check2 == 1
    % Rotate the orientations
    rTex = fabOr.*o;
    warning on
    
else
    fprintf('\n\n!!!  ERROR !!!\nOne of the input variables is *not* an orientation\nNote: If you are trying to rotate an ODF, please use rotODFGeo function\n')
    rTex = [];
end