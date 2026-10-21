function test_geographic_convention_proof
% Check a public, synthetic field row through all eight cardinal layouts.
here = fileparts(mfilename('fullpath'));
csvPath = fullfile(here,'field_fabric_synthetic.csv');
field = readtable(csvPath,'TextType','string');
assert(height(field) == 1);
assert(strcmp(field.name(1),"SYNTHETIC DOWN-DIP"));
assert(isequal([field.strike,field.dip,field.trend,field.plunge], ...
    [30,60,120,60]));

% The no-argument defaults must be independently runnable after a clone.
cases = geographic_convention_proof([],[],false);
explicitCases = geographic_convention_proof(csvPath,'SYNTHETIC DOWN-DIP',false);
try
    geographic_convention_proof(csvPath,'SYNTHETIC DOWN-DIP',true);
    error('test_geographic_convention_proof:MissingGuard', ...
        'A custom CSV unexpectedly used the bundled synthetic output path.');
catch ME
    assert(strcmp(ME.identifier,'geographic_convention_proof:OutputRequired'), ...
        'A custom CSV must require an explicit plot output path.');
end
assert(numel(cases) == 8);
assert(numel(explicitCases) == 8);
for k = 1:numel(cases)
    c = cases(k);
    assert(strcmp(c.sampleName,"SYNTHETIC DOWN-DIP"));
    assert(isequal([c.strike,c.dip,c.trend,c.plunge],[30,60,120,60]));
    assert(c.planeResidual < 1e-10);
    assert(c.measuredLineAngleDeg < 1e-5);
    assert(c.mappedLineAngleDeg < 1e-5);
    assert(c.mappedPoleAngleDeg < 1e-5);
    assert(max(abs(c.geoLine - cases(1).geoLine)) < 1e-10);
    assert(max(abs(c.geoPole - cases(1).geoPole)) < 1e-10);
    assert(max(abs(c.geoLine - explicitCases(k).geoLine)) < 1e-10);
    assert(max(abs(c.geoPole - explicitCases(k).geoPole)) < 1e-10);
end
fprintf('Eight synthetic geographic convention cases passed.\n');
end
