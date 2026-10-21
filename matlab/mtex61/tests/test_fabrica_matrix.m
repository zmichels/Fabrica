function report = test_fabrica_matrix
%TEST_FABRICA_MATRIX Deterministic spherical-grid audit of geographic fabric.
% The independent oracle is numerical east/north/down geometry. It does not
% use Fabrica to construct expected field directions. The grid crosses
% specimen line/pole pairs, valid strike/dip/pitch field cases, and eight
% cardinal MTEX plotting conventions. Run after startup_mtex and addpath.

started = tic;

% Equal-area-like coverage of the specimen sphere, plus exact cardinal axes.
nFibonacci = 24;
specimenLines = zeros(nFibonacci+6,3);
goldenAngle = pi*(3-sqrt(5));
for i = 0:nFibonacci-1
    z = 1 - 2*(i+0.5)/nFibonacci;
    azimuth = goldenAngle*i;
    radius = sqrt(1-z*z);
    specimenLines(i+1,:) = [radius*cos(azimuth),radius*sin(azimuth),z];
end
specimenLines(nFibonacci+1:end,:) = [eye(3);-eye(3)];
poleAngles = [0,45,90,135];

% Right-hand-rule strike, dip, and pitch from strike toward down-dip.
% A vertical foliation's trend cannot determine a nonhorizontal plunge,
% so the main trend-only matrix keeps its two horizontal strike directions.
strikeValues = [0,37,89,137,179,225,271,322];
dipValues = [0,0.01,1,19,45,71,89,89.99,90];
pitchValues = [0,0.1,15,45,90,135,165,179.9,180];
[strike,dip,pitch] = ndgrid(strikeValues,dipValues,pitchValues);
keep = ~(dip == 90 & pitch ~= 0 & pitch ~= 180);
strike = strike(keep);
dip = dip(keep);
pitch = pitch(keep);
nField = numel(strike);

% Numeric ENDown oracle: strike S, horizontal dip H, down-dip Q,
% downward plane pole P, and the in-plane lineation L.
S = [sind(strike),cosd(strike),zeros(nField,1)];
H = [cosd(strike),-sind(strike),zeros(nField,1)];
Q = cosd(dip).*H + [zeros(nField,2),sind(dip)];
expectedPole = -sind(dip).*H + [zeros(nField,2),cosd(dip)];
expectedLine = cosd(pitch).*S + sind(pitch).*Q;
expectedTrend = mod(atan2d(expectedLine(:,1),expectedLine(:,2)),360);
expectedPlunge = asind(expectedLine(:,3));
expectedLineXY = expectedLine(:,1:2)./sqrt(1+expectedLine(:,3));
expectedPoleXY = expectedPole(:,1:2)./sqrt(1+expectedPole(:,3));
assert(max(abs(sum(expectedLine.*expectedPole,2))) < 1e-12);

outAxes = {zvector,-zvector};
outLabels = {'+Z','-Z'};
eastAxes = {xvector,yvector,-xvector,-yvector};
eastLabels = {'+X','+Y','-X','-Y'};
names = {'strikeComponent','lineComponent','poleComponent', ...
    'mappedLineComponent','mappedPoleComponent','mappedThirdComponent', ...
    'planeDot','trendDeg','plungeDeg','equalAreaXY'};
worst = struct();
for i = 1:numel(names)
    worst.(names{i}) = struct('value',-Inf,'fieldIndex',0, ...
        'specimenIndex',0,'poleAngle',0,'screenOut','', ...
        'screenEast','','strike',0,'dip',0,'pitch',0, ...
        'trend',0,'plunge',0);
end
evaluations = 0;
specimenPairs = size(specimenLines,1)*numel(poleAngles);
totalEvaluations = specimenPairs*nField*numel(outAxes)*numel(eastAxes);
directionResiduals = zeros(totalEvaluations,1);
projectionResiduals = zeros(totalEvaluations,1);

for iz = 1:numel(outAxes)
    for ie = 1:numel(eastAxes)
        pc = plottingConvention(outAxes{iz},eastAxes{ie});
        for il = 1:size(specimenLines,1)
            lineNumeric = specimenLines(il,:);
            if abs(lineNumeric(3)) < 0.9
                reference = [0,0,1];
            else
                reference = [0,1,0];
            end
            tangentA = cross(lineNumeric,reference);
            tangentA = tangentA/norm(tangentA);
            tangentB = cross(lineNumeric,tangentA);
            specimenLine = vector3d(lineNumeric(1),lineNumeric(2),lineNumeric(3));

            for ia = 1:numel(poleAngles)
                poleNumeric = cosd(poleAngles(ia))*tangentA + ...
                    sind(poleAngles(ia))*tangentB;
                specimenPole = vector3d(poleNumeric(1),poleNumeric(2),poleNumeric(3));
                assert(abs(dot(specimenLine,specimenPole)) < 1e-12);

                [~,fabRot,strikeV,poleV,lineV] = geoFab( ...
                    strike,dip,expectedTrend,specimenLine,specimenPole,pc);
                gotStrike = geographicComponents(strikeV,pc);
                gotLine = geographicComponents(lineV,pc);
                gotPole = geographicComponents(poleV,pc);
                mappedLine = geographicComponents(fabRot*specimenLine,pc);
                mappedPole = geographicComponents(fabRot*specimenPole,pc);
                sourceThird = cross(specimenLine,specimenPole);
                mappedThird = fabRot*sourceThird;
                targetThird = cross(lineV,poleV);
                thirdError = sqrt((mappedThird.x-targetThird.x).^2 + ...
                    (mappedThird.y-targetThird.y).^2 + ...
                    (mappedThird.z-targetThird.z).^2);
                [gotTrend,gotPlunge] = V2TP(lineV,pc);
                gotLineXY = gotLine(:,1:2)./sqrt(1+gotLine(:,3));
                gotPoleXY = gotPole(:,1:2)./sqrt(1+gotPole(:,3));
                mappedLineXY = mappedLine(:,1:2)./sqrt(1+mappedLine(:,3));
                mappedPoleXY = mappedPole(:,1:2)./sqrt(1+mappedPole(:,3));

                metrics = struct( ...
                    'strikeComponent',vecnorm(gotStrike-S,2,2), ...
                    'lineComponent',vecnorm(gotLine-expectedLine,2,2), ...
                    'poleComponent',vecnorm(gotPole-expectedPole,2,2), ...
                    'mappedLineComponent',vecnorm(mappedLine-expectedLine,2,2), ...
                    'mappedPoleComponent',vecnorm(mappedPole-expectedPole,2,2), ...
                    'mappedThirdComponent',thirdError(:), ...
                    'planeDot',abs(dot(poleV,lineV)), ...
                    'trendDeg',abs(mod(gotTrend-expectedTrend+180,360)-180), ...
                    'plungeDeg',abs(gotPlunge-expectedPlunge), ...
                    'equalAreaXY',max(vecnorm(gotLineXY-expectedLineXY,2,2), ...
                        vecnorm(gotPoleXY-expectedPoleXY,2,2)));
                indices = evaluations+(1:nField);
                directionResiduals(indices) = max( ...
                    metrics.mappedLineComponent,metrics.mappedPoleComponent);
                projectionResiduals(indices) = max( ...
                    vecnorm(mappedLineXY-expectedLineXY,2,2), ...
                    vecnorm(mappedPoleXY-expectedPoleXY,2,2));
                for im = 1:numel(names)
                    name = names{im};
                    values = metrics.(name);
                    [value,index] = max(values(:));
                    if value > worst.(name).value
                        worst.(name) = struct('value',value, ...
                            'fieldIndex',index,'specimenIndex',il, ...
                            'poleAngle',poleAngles(ia), ...
                            'screenOut',outLabels{iz},'screenEast',eastLabels{ie}, ...
                            'strike',strike(index),'dip',dip(index), ...
                            'pitch',pitch(index),'trend',expectedTrend(index), ...
                            'plunge',expectedPlunge(index));
                    end
                end
                evaluations = evaluations+nField;
            end
        end
    end
end

% A trend is undefined for a vertical lineation. The pitch-based entry
% point retains this exact case, which the trend-only matrix omits.
verticalPitchCases = 0;
for iz = 1:numel(outAxes)
    for ie = 1:numel(eastAxes)
        pc = plottingConvention(outAxes{iz},eastAxes{ie});
        [~,~,~,verticalPole,verticalLine] = SDP2or(37,90,90,pc);
        [verticalTrend,verticalPlunge] = V2TP(verticalLine,pc);
        assert(isnan(verticalTrend) && abs(verticalPlunge-90) < 1e-8);
        assert(abs(dot(verticalPole,verticalLine)) < 1e-10);
        verticalPitchCases = verticalPitchCases+1;
    end
end

% Deliberately read an asymmetric, correctly constructed field vector with
% the wrong screen convention. A wrong pole side must fail plane membership.
goodPc = plottingConvention(zvector,xvector);
badPc = plottingConvention(zvector,-xvector);
controlStrike = 37;
controlDip = 48;
controlPitch = 31;
controlS = [sind(controlStrike),cosd(controlStrike),0];
controlH = [cosd(controlStrike),-sind(controlStrike),0];
controlQ = cosd(controlDip)*controlH + [0,0,sind(controlDip)];
controlExpected = cosd(controlPitch)*controlS + sind(controlPitch)*controlQ;
controlTrend = mod(atan2d(controlExpected(1),controlExpected(2)),360);
controlPlunge = asind(controlExpected(3));
controlLine = TP2vec(controlTrend,controlPlunge,goodPc);
goodGeo = geographicComponents(controlLine,goodPc);
badGeo = geographicComponents(controlLine,badPc);
wrongConventionGap = norm(goodGeo-badGeo);
assert(wrongConventionGap > 1);
wrongPole = TP2vec(controlStrike+90,90-controlDip,goodPc);
wrongPolePlaneDot = abs(dot(wrongPole,controlLine));
assert(wrongPolePlaneDot > 0.1);

tolerance = struct('components',1e-9,'planeDot',1e-10, ...
    'anglesDeg',1e-6,'equalAreaXY',1e-9);
for im = 1:numel(names)
    name = names{im};
    if strcmp(name,'planeDot')
        limit = tolerance.planeDot;
    elseif strcmp(name,'trendDeg') || strcmp(name,'plungeDeg')
        limit = tolerance.anglesDeg;
    elseif strcmp(name,'equalAreaXY')
        limit = tolerance.equalAreaXY;
    else
        limit = tolerance.components;
    end
    assert(worst.(name).value < limit, ...
        '%s exceeded tolerance: %.12g >= %.12g', ...
        name,worst.(name).value,limit);
end
assert(evaluations == totalEvaluations);
assert(max(directionResiduals) < tolerance.components);
assert(max(projectionResiduals) < tolerance.equalAreaXY);
histogramEdges = -17:0.5:-8;
residualDistributions = struct( ...
    'directionUnitVector',residualDistribution(directionResiduals,histogramEdges), ...
    'equalAreaUnitHorizon',residualDistribution(projectionResiduals,histogramEdges));

report = struct('specimenLineDirections',size(specimenLines,1), ...
    'specimenPoleAngles',poleAngles,'specimenPairs',specimenPairs, ...
    'fieldCases',nField,'plottingConventions',8, ...
    'verticalPitchCases',verticalPitchCases, ...
    'evaluations',evaluations,'elapsedSeconds',toc(started), ...
    'tolerance',tolerance,'worst',worst, ...
    'residualDistributions',residualDistributions, ...
    'residualMetricDefinitions',struct( ...
        'directionUnitVector', ...
        'For each mapping, max Euclidean ENDown unit-vector error of rotated specimen line and pole against the independent field oracle.', ...
        'equalAreaUnitHorizon', ...
        'For each mapping, max Euclidean equal-area XY error of rotated specimen line and pole against the independent field oracle. XY=(East,North)/sqrt(1+Down); horizon radius is 1.', ...
        'histogramFloor', ...
        'Errors below 1e-17, including exact zeros, are placed in the first log10 histogram bin.'), ...
    'negativeControlWrongConventionGap',wrongConventionGap, ...
    'negativeControlWrongPolePlaneDot',wrongPolePlaneDot);
here = fileparts(mfilename('fullpath'));
reportPath = fullfile(here,'..','geographic_matrix_report_mtex61.json');
fid = fopen(reportPath,'w');
assert(fid ~= -1,'Cannot write matrix report.');
cleaner = onCleanup(@() fclose(fid));
fwrite(fid,jsonencode(report));
fprintf('Fabrica matrix: %d specimen pairs x %d field cases x 8 conventions = %d mappings.\n', ...
    specimenPairs,nField,evaluations);
fprintf('Worst line %.3g, pole %.3g, rotation %.3g, plane dot %.3g, XY %.3g; %.2f s.\n', ...
    worst.lineComponent.value,worst.poleComponent.value, ...
    max(worst.mappedLineComponent.value,worst.mappedPoleComponent.value), ...
    worst.planeDot.value,worst.equalAreaXY.value,report.elapsedSeconds);
end

function components = geographicComponents(v,pc)
components = [dot(v,pc.east),dot(v,pc.north),dot(v,pc.intoScreen)];
components = reshape(components,[],3);
end

function distribution = residualDistribution(values,edges)
values = sort(values(:));
n = numel(values);
probabilities = [0.5,0.9,0.99,0.999];
positions = 1+(n-1)*probabilities;
low = floor(positions);
high = ceil(positions);
fraction = positions-low;
quantileValues = (1-fraction).*values(low).' + fraction.*values(high).';
quantiles = struct('p50',quantileValues(1),'p90',quantileValues(2), ...
    'p99',quantileValues(3),'p99_9',quantileValues(4), ...
    'max',values(end));
logValues = log10(max(values,1e-17));
counts = histcounts(logValues,edges);
assert(sum(counts) == n,'Residual histogram does not cover all mappings.');
distribution = struct('count',n,'quantiles',quantiles, ...
    'log10BinEdges',edges,'log10BinCounts',counts);
end
