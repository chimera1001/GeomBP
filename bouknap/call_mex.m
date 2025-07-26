clc
clear all


w = int64([3,4,5,6,7]);
p = int64([4,6,5,6,7]);
d = int32([0,1,0,0,0]);
c = int64(20);
n = int32(5);

[x, z] = bouknap(w, p, d, c, n);