%% Synthetic teaching pairs in example.csv: dip direction, dip, right-hand rake
import_example;

% This convention places geographic east at +X and north at +Y.
pc = plottingConvention(zvector,xvector);
[fabOr, strikes, poles, lins] = FandL2Or(dd,d,r,pc); %#ok<ASGLU>

%% Density and mode of the fabric orientations
fabODF = calcDensity(fabOr,'halfwidth',5*degree);
[~,fabMode] = max(fabODF);
cs = fabOr.CS;
h = [Miller(1,0,0,'direction',cs),Miller(0,0,1,'direction',cs)];

%% Lower-hemisphere pole figures in the same geographic screen frame
figure;
plotPDF(fabODF,h,'antipodal','lower','smooth','colorRange','equal',pc);
hold on;
plotPDF(fabMode,h,'antipodal','lower','Marker','s','MarkerSize',12, ...
    'MarkerFaceColor','k','add2all');
annotate([pc.north,pc.east,pc.south,pc.west], ...
    'label',{'N','E','S','W'},'BackgroundColor','w');
