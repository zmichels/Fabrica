function test_randFab
%TEST_RANDFAB Check geographic random-pair output against MTEX fabric vectors.

previousRng = rng;
restoreRng = onCleanup(@()rng(previousRng)); %#ok<NASGU>
rng(261020);
[s,d,dd,t,p,r] = randFab(128);
assert(all(abs(mod(dd-s-90+180,360)-180) < 1e-12));
assert(any(r > 90) && any(r < 90));
assert(all(p >= 0 & p <= 90));

pc = plottingConvention(zvector,xvector);
[~,~,~,~,line] = SDP2or(s,d,r,pc);
[gotTrend,gotPlunge] = V2TP(line,pc);
assert(all(abs(mod(gotTrend-t+180,360)-180) < 1e-7));
assert(all(abs(gotPlunge-p) < 1e-7));

% A >90-degree rake must point into the opposite-strike half of the plane.
index = find(r > 120,1);
assert(~isempty(index));
assert(cosd(t(index)-s(index)) < 0);
fprintf('randFab right-hand-rake quadrants passed for %d pairs.\n',numel(s));
end
