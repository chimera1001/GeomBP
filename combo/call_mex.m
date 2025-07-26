clc
clear all

w = int64([3,4,5,6,7]);
p = int64([4,6,5,6,7]);
c = int64(12);
lb = int64(0);
ub = int64(0);
def = int32(1);
relx = int32(0);
n = int32(5);

[sol,z] = combo(w,p,c,lb,ub,def,relx,n);