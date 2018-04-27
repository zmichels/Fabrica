function [s,d,dd,t,p,r] = randFab(num)

%%

% random fabric elements
s = rand(num,1)*360;
dd = s+90;
dd(dd>360) = dd(dd>360)-360;
d = rand(num,1)*90;
r = rand(num,1)*180;

% trend and plunge from rake (right hand rule, 0-180) 
% trend
t = s + atand((tand(r).*cosd(d)));


% plunge
p = asind(sind(d).*sind(r));

% correction for negatives
t(p<0) = t(p<0)+180;
t(t<0) = t(t<0)+360;
t(t>360) = t(t>360)-360;

p(p<0) = -p(p<0);


