function test_fabrica_mtex61
% Run after startup_mtex and addpath on Fabrica and this directory.
outs = {zvector,-zvector};
easts = {xvector,yvector,-xvector,-yvector};
trends = [40,130,220,310];
plunge = 30;
for iz = 1:numel(outs)
  for ie = 1:numel(easts)
    pc = plottingConvention(outs{iz},easts{ie});
    for trend = trends
      v = TP2vec(trend,plunge,pc);
      expected = cosd(plunge) * (sind(trend) * pc.east + ...
        cosd(trend) * pc.north) + sind(plunge) * pc.intoScreen;
      assert(angle(v,expected,'noSymmetry') < 1e-6, ...
        'Wrong trend/plunge for screen convention %d,%d, trend %g', ...
        iz,ie,trend);
      [gotTrend,gotPlunge] = V2TP(v,pc);
      assert(abs(mod(gotTrend-trend+180,360)-180) < 1e-6);
      assert(abs(gotPlunge-plunge) < 1e-6);
    end
    [verticalTrend,verticalPlunge] = V2TP(TP2vec(123,90,pc),pc);
    assert(isnan(verticalTrend));
    assert(abs(verticalPlunge-90) < 1e-6);
    [fabOr,fabRot,~,poleV,linV] = SDT2or(15,56,42,pc);
    if dot(pc.outOfScreen,zvector) > 0
      linDir = xvector;
      poleDir = -yvector;
    else
      linDir = -xvector;
      poleDir = yvector;
    end
    assert(angle(fabRot*linDir,linV,'noSymmetry') < 1e-6);
    assert(angle(fabRot*poleDir,poleV,'noSymmetry') < 1e-6);
    assert(abs(dot(linV,poleV)) < 1e-6);
    assert(fabOr.CS == crystalSymmetry('mmm'));
    [gotTrend,~] = V2TP(linV,pc);
    assert(abs(mod(gotTrend-42+180,360)-180) < 1e-6);

    [~,pitchRot,~,pitchPole,pitchLin] = SDP2or(15,56,25,pc);
    assert(angle(pitchRot*linDir,pitchLin,'noSymmetry') < 1e-6);
    assert(angle(pitchRot*poleDir,pitchPole,'noSymmetry') < 1e-6);

    [rakeOr,strikeV,rakePole,rakeLin] = FandL2Or(90,45,90,pc);
    assert(angle(strikeV,TP2vec(0,0,pc),'noSymmetry') < 1e-6);
    assert(angle(rakePole,TP2vec(270,45,pc),'noSymmetry') < 1e-6);
    assert(angle(rakeLin,TP2vec(90,45,pc),'noSymmetry') < 1e-6);
    assert(angle(rakeOr*zvector,rakePole,'noSymmetry') < 1e-6);
    assert(angle(rakeOr*xvector,rakeLin,'noSymmetry') < 1e-6);

    sample = makeGeoSample(15,56,42,pc);
    assert(angle(sample.spec.lin,linDir,'noSymmetry') < 1e-6);
    assert(angle(sample.spec.pole,poleDir,'noSymmetry') < 1e-6);
    [~,linePlunge] = V2TP(sample.geo.lin,pc);
    assert(abs(sample.geo.plunge-linePlunge) < 1e-6);
  end
end

pc = plottingConvention(zvector,xvector);
% Vertical foliation admits strike-parallel lineation without a unique
% plunge; an off-strike trend is geometrically incompatible.
[~,verticalRot,~,verticalPole,verticalLin] = SDT2or(0,90,0,pc);
assert(angle(verticalRot*xvector,verticalLin,'noSymmetry') < 1e-6);
assert(abs(dot(verticalLin,verticalPole)) < 1e-6);
didError = false;
try
    SDT2or(0,90,90,pc);
catch ME
    didError = strcmp(ME.identifier,'geoFab:Trend');
end
assert(didError);

% The included synthetic CSV supplies dip-direction, dip, and rake arrays.
here = fileparts(mfilename('fullpath'));
pairs = readtable(fullfile(here,'..','example.csv'));
[manyOr,~,manyPoles,manyLins] = FandL2Or(pairs.dipDir,pairs.dip,pairs.rhrake,pc);
assert(numel(manyOr) == height(pairs));
assert(all(abs(dot(manyPoles,manyLins)) < 1e-6));
run(fullfile(here,'..','import_example.m'));
assert(numel(dd) == height(pairs));

pc = plottingConvention(zvector,xvector);
[~,fabRot] = SDT2or(15,56,42,pc);
cs = crystalSymmetry('mmm');
ori = orientation.byEuler(10*degree,20*degree,30*degree,cs);
h = Miller(1,0,0,'direction',cs);
rotOri = rotTex2Geo(ori,fabRot);
assert(angle(rotOri*h,fabRot*(ori*h),'noSymmetry') < 1e-6);

pos = vector3d([0;1;0;1],[0;0;1;1],zeros(4,1));
ebsd = EBSD(pos,rotation.id(4,1),ones(4,1),{cs},struct());
rotEbsd = rotEBSD2Geo(ebsd,fabRot);
assert(all(abs(rotEbsd.x-ebsd.x) < 1e-10));
assert(all(abs(rotEbsd.y-ebsd.y) < 1e-10));
assert(angle(rotEbsd('indexed').orientations(1)*h, ...
    fabRot*(ebsd('indexed').orientations(1)*h),'noSymmetry') < 1e-6);

sample = makeGeoSample(15,56,42,ebsd,pc);
assert(angle(sample.geo.rotebsd{1}('indexed').orientations(1)*h, ...
    fabRot*(ebsd('indexed').orientations(1)*h),'noSymmetry') < 1e-6);
updipSample = makeGeoSample(0,60,270,pc);
assert(updipSample.geo.trend == 270);
assert(abs(updipSample.geo.lineTrend-90) < 1e-6);

odf = unimodalODF(ori);
rotOdf = rotODFGeo(odf,fabRot);
assert(isa(rotOdf,'SO3Fun'));
assert(abs(eval(rotOdf,rotOri)-eval(odf,ori)) < 1e-5);
fprintf('Fabrica MTEX 6.1 direction and fabric cases passed.\n');
end
