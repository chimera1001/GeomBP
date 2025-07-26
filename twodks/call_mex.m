clc
clear all

w = int32([3,4,5,6,7]);
p = int32([4,6,5,6,7]);
cw = 12;
cu = 14;


[sol, z] = twodks(p, w, p, cw, cu);