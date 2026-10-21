function runFabricaTests
%RUNFABRICATESTS Verify the public MATLAB workflow with MTEX 6.1.1 or 7.2.
% Run startup_mtex first, then add the Fabrica directory to the MATLAB path.

here = fileparts(mfilename('fullpath'));
addpath(here,fullfile(here,'tests'));
if exist('plottingConvention','class') ~= 8
    error('runFabricaTests:MTEX', ...
        'MTEX is not initialized. Run startup_mtex before this suite.');
end

test_fabrica_mtex61;
test_fabrica_wrappers;
test_randFab;
test_fabrica_plot_frames;
test_geographic_convention_proof;
test_fabrica_matrix;
fprintf('All Fabrica MATLAB tests passed under %s.\n',getMTEXpref('version'));
end
