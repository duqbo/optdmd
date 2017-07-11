
% build mex binaries

fortheaderdir = fullfile(matlabroot,'extern','examples','refbook');
fortfile = fullfile(matlabroot,'extern','examples','refbook','fort.c');

mex('-v','-outdir','./bin','-largeArrayDims',['-I' fortheaderdir],'./src/xgeqp3_m.c',fortfile,'-lmwlapack')
mex('-v','-outdir','./bin','-largeArrayDims',['-I' fortheaderdir],'./src/xormqr_m.c',fortfile,'-lmwlapack')
