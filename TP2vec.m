function [v] = TP2vec(trend,plunge)
v = vector3d('polar',(90+plunge)*degree,(90-trend)*degree);
