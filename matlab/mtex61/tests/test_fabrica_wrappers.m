function test_fabrica_wrappers
%TEST_FABRICA_WRAPPERS Cross-check the three public fabric entry points.
% Run after startup_mtex and addpath on Fabrica and this directory.

fields = [ ... % RHR strike, dip, rake from strike toward down dip
      0,  0,   0;
     15, 56,  25;
     89, 89.99, 135;
    225, 45,  90;
    322, 71, 179.9];
outs = {zvector,-zvector};
easts = {xvector,yvector,-xvector,-yvector};

for iz = 1:numel(outs)
    for ie = 1:numel(easts)
        pc = plottingConvention(outs{iz},easts{ie});
        for i = 1:size(fields,1)
            s = fields(i,1); d = fields(i,2); r = fields(i,3);
            [~,rotP,strikeP,poleP,lineP] = SDP2or(s,d,r,pc);
            [trend,~] = V2TP(lineP,pc);
            [~,rotT,strikeT,poleT,lineT] = SDT2or(s,d,trend,pc);
            [orF,strikeF,poleF,lineF] = FandL2Or(mod(s+90,360),d,r,pc);

            assert(vectorError(strikeP,strikeT) < 1e-9);
            assert(vectorError(poleP,poleT) < 1e-9);
            assert(vectorError(lineP,lineT) < 1e-9);
            assert(vectorError(strikeP,strikeF) < 1e-9);
            assert(vectorError(poleP,poleF) < 1e-9);
            assert(vectorError(lineP,lineF) < 1e-9);

            if dot(pc.outOfScreen,zvector) > 0
                specimenLine = xvector;
                specimenPole = -yvector;
            else
                specimenLine = -xvector;
                specimenPole = yvector;
            end
            assert(vectorError(rotP*specimenLine,lineP) < 1e-9);
            assert(vectorError(rotP*specimenPole,poleP) < 1e-9);
            assert(vectorError(rotT*specimenLine,lineT) < 1e-9);
            assert(vectorError(rotT*specimenPole,poleT) < 1e-9);
            assert(vectorError(orF*xvector,lineF) < 1e-9);
            assert(vectorError(orF*zvector,poleF) < 1e-9);
            assert(abs(dot(lineP,poleP)) < 1e-10);
        end
    end
end

function error = vectorError(a,b)
% Cartesian differences stay accurate when inverse-cosine angles round off.
error = norm([a.x-b.x,a.y-b.y,a.z-b.z]);
end
fprintf('Fabric wrapper equivalence passed for %d field/convention cases.\n', ...
    size(fields,1)*numel(outs)*numel(easts));
end
