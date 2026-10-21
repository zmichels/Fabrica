function test_fabrica_plot_frames
% Verify the specimen axes used by the shared fabric plotting helper.
outs = {zvector,-zvector};
easts = {xvector,yvector,-xvector,-yvector};
for i = 1:numel(outs)
    for j = 1:numel(easts)
        pc = plottingConvention(outs{i},easts{j});
        [fabOr,~,~,poleV,lineV] = SDT2or(15,56,42,pc);
        if dot(pc.outOfScreen,zvector) > 0
            sourceLine = xvector;
            sourcePole = -yvector;
        else
            sourceLine = -xvector;
            sourcePole = yvector;
        end
        [plotLine,plotPole] = resolveFabPlotVectors( ...
            fabOr,sourceLine,sourcePole);
        assert(angle(plotLine,lineV,'noSymmetry') < 1e-6);
        assert(angle(plotPole,poleV,'noSymmetry') < 1e-6);

        [legacyOr,~,legacyPole,legacyLine] = FandL2Or(90,45,90,pc);
        [plotLine,plotPole] = resolveFabPlotVectors(legacyOr);
        assert(angle(plotLine,legacyLine,'noSymmetry') < 1e-6);
        assert(angle(plotPole,legacyPole,'noSymmetry') < 1e-6);
    end
end
try
    compareInputs = nargin('plot_compare');
catch
    compareInputs = NaN;
end
assert(compareInputs == -6, ...
    'plot_compare must be callable with explicit data and convention inputs.');
fprintf('Fabrica plot source-frame cases passed.\n');
end
