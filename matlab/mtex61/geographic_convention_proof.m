function cases = geographic_convention_proof(csvPath,sampleName,makePlot,outputPath)
%GEOGRAPHIC_CONVENTION_PROOF Check one field fabric in eight screen layouts.
% The default CSV contains a public synthetic down-dip example. An explicit
% CSV may contain measured data. Both use right-hand-rule strike and downward
% dip, plus lineation trend and plunge (degrees). Each case uses a different
% cardinal plotting convention and orthogonal specimen line/pole axes.
% With makePlot=true, MTEX draws eight geographic lower-hemisphere panels.

customCsv = nargin >= 1 && ~isempty(csvPath);
if ~customCsv
    csvPath = fullfile(fileparts(mfilename('fullpath')), ...
        'tests','field_fabric_synthetic.csv');
end
if nargin < 2 || isempty(sampleName)
    sampleName = 'SYNTHETIC DOWN-DIP';
end
if nargin < 3 || isempty(makePlot)
    makePlot = true;
end
if nargin < 4 || isempty(outputPath)
    if customCsv && makePlot
        error('geographic_convention_proof:OutputRequired', ...
            'Supply outputPath for a custom CSV so the bundled synthetic figure is preserved.');
    end
    outputPath = fullfile(fileparts(mfilename('fullpath')), ...
        'geographic_convention_proof_mtex61.png');
end

field = readtable(csvPath,'TextType','string');
row = field(strcmp(strtrim(string(field.name)),string(sampleName)),:);
assert(height(row) == 1,'Expected exactly one field-data row for %s.',sampleName);
strike = row.strike;
dip = row.dip;
trend = row.trend;
plunge = row.plunge;
dipDirection = strtrim(string(row.dipDir));
assert(isscalar(strike) && isscalar(dip) && isscalar(trend) && isscalar(plunge));
assert(all(isfinite([strike,dip,trend,plunge])));

screenOut = {zvector,-zvector};
outLabels = {'+Z','-Z'};
screenEast = {xvector,yvector,-xvector,-yvector};
eastLabels = {'+X','+Y','-X','-Y'};
% The eight pairs illustrate different perpendicular specimen-axis
% assignments. The first uses +X lineation and +Y foliation pole.
specLine = {xvector,yvector,-xvector,-yvector, ...
    xvector,yvector,-xvector,-yvector};
specPole = {yvector,-xvector,-yvector,xvector, ...
    -yvector,xvector,yvector,-xvector};
lineLabels = {'+X','+Y','-X','-Y','+X','+Y','-X','-Y'};
poleLabels = {'+Y','-X','-Y','+X','-Y','+X','+Y','-X'};

cases = struct([]);
for iz = 1:2
    for ie = 1:4
        k = (iz-1)*4 + ie;
        pc = plottingConvention(screenOut{iz},screenEast{ie});
        measuredPole = TP2vec(strike-90,90-dip,pc);
        measuredLine = TP2vec(trend,plunge,pc);
        [~,fabRot,~,poleV,linV] = geoFab( ...
            strike,dip,trend,specLine{k},specPole{k},pc);
        geoLine = [dot(linV,pc.east),dot(linV,pc.north), ...
            dot(linV,pc.intoScreen)];
        geoPole = [dot(poleV,pc.east),dot(poleV,pc.north), ...
            dot(poleV,pc.intoScreen)];
        oneCase = struct( ...
            'sampleName',string(sampleName),'strike',strike,'dip',dip, ...
            'trend',trend,'plunge',plunge,'pc',pc, ...
            'outLabel',outLabels{iz},'eastLabel',eastLabels{ie}, ...
            'specimenLineLabel',lineLabels{k}, ...
            'specimenPoleLabel',poleLabels{k}, ...
            'pole',poleV,'lineation',linV,'measuredLineation',measuredLine, ...
            'rotation',fabRot,'geoLine',geoLine,'geoPole',geoPole, ...
            'planeResidual',abs(dot(measuredLine,measuredPole)), ...
            'measuredLineAngleDeg',angle(linV,measuredLine,'noSymmetry')/degree, ...
            'mappedLineAngleDeg',angle(fabRot*specLine{k},measuredLine,'noSymmetry')/degree, ...
            'mappedPoleAngleDeg',angle(fabRot*specPole{k},measuredPole,'noSymmetry')/degree);
        if k == 1
            cases = oneCase;
        else
            cases(k) = oneCase;
        end
        assert(cases(k).planeResidual < 1e-10, ...
            'Input lineation does not lie in the input foliation.');
        assert(cases(k).measuredLineAngleDeg < 1e-5, ...
            'Computed lineation differs from the input trend/plunge.');
    end
end

if ~makePlot
    return
end

newMtexFigure('layout',[3,4],'figSize','huge');
for k = 1:numel(cases)
    if k > 1
        if k <= 4
            nextAxis(1,k);
        else
            nextAxis(3,k-4);
        end
    end
    c = cases(k);
    plot(c.pole,'antipodal','lower','plane',c.pc, ...
        'LineColor',[0.15 0.35 0.65],'LineWidth',1.8);
    hold on
    % MTEX annotate defaults to add2all; bind each symbol to this panel.
    annotate(c.pole,'Marker','s','MarkerSize',11, ...
        'MarkerEdgeColor',[0.05 0.25 0.55], ...
        'MarkerFaceColor',[0.05 0.25 0.55],'parent',gca);
    annotate(c.measuredLineation,'Marker','o','MarkerSize',10, ...
        'MarkerEdgeColor',[0.75 0.18 0.12], ...
        'MarkerFaceColor',[0.75 0.18 0.12],'parent',gca);
    annotate([c.pc.north,c.pc.east,c.pc.south,c.pc.west], ...
        'label',{'N','E','S','W'},'Marker','none', ...
        'BackgroundColor','w','parent',gca);
    hold off
    % Keep convention and specimen labels outside the stereonet so they
    % cannot cover the foliation trace or its line/pole symbols.
    column = mod(k-1,4);
    if k <= 4
        labelY = 0.64;
    else
        labelY = 0.34;
    end
    annotation(gcf,'textbox',[0.015+0.25*column,labelY,0.22,0.05], ...
        'String',sprintf('%s out, %s east\nline %s, pole %s', ...
        c.outLabel,c.eastLabel,c.specimenLineLabel,c.specimenPoleLabel), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',10,'EdgeColor','none');
end
annotation(gcf,'textbox',[0.08,0.48,0.84,0.055], ...
    'String',sprintf('%s: RHR strike %.0f deg, dip %.0f deg (%s); lineation trend %.0f deg, plunge %.0f deg', ...
    sampleName,strike,dip,char(dipDirection),trend,plunge), ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',15,'EdgeColor','none');
annotation(gcf,'textbox',[0.08,0.45,0.84,0.03], ...
    'String','Blue = foliation great circle   |   square = foliation pole   |   circle = input lineation', ...
    'HorizontalAlignment','center','FontSize',11,'EdgeColor','none');
exportgraphics(gcf,outputPath,'Resolution',180);
end
