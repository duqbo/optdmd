
%% Set up for optdmd
%% open a MatLab session and run this script
%% in this folder first,
%% before running the scripts in the examples 
%% folder

% add appropriate folders to path

addpath('./src');
addpath('./examples');

% make sure that the lapack mex file has been
% compiled. respond 'y' to install

lapack('dgesvd');

