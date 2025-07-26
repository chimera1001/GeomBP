clc
clear all

mex COPTIMFLAGS="-O2" twodks.c twodks.o -lgfortran

% mex COPTIMFLAGS="-O2" combo.c