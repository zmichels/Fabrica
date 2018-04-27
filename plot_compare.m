%% plots to compare and visually check that rotation is working

f1 = figure;
hold on
plotPDF(o,h,'antipodal','lower','smooth','colorrange','equal')
annotate(xvector,'antipodal','lower','Marker','d','MarkerZise',12,'MarkerFaceColor','k')
annotate(yvector,'antipodal','lower','plane')
hold off


f2 = figure;
hold on
plotPDF(rTex,h,'antipodal','lower','smooth','colorrange','equal')
annotate(fabOr*xvector,'antipodal','lower','Marker','d','MarkerZise',12,'MarkerFaceColor','k')
annotate(fabOr*zvector,'antipodal','lower','plane')