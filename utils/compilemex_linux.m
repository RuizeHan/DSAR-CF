% Run this file to build the needed mex-files on windows

% Build merResize
mex -lopencv_core -lopencv_imgproc -L./ -I./ mexResize.cpp MxArray.cpp

% Build setnonzeros from the lightspeed matlab toolbox
mex setnonzeros.c

% Build gradientMex from Piotrs toolbox
mex gradientMex.cpp -I./

% Build mtimesx
mex -DDEFINEUNIX -largeArrayDims -L/usr/local/MATLAB/R2015a/bin/glnxa64/ -lmwblas mtimesx.c