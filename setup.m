
%% Set up for optdmd
%% open a MatLab session and run this script
%% in this folder first,
%% before running the scripts in the examples 
%% folder

% create bin folder if it doesn't exist

if ~exist('./bin', 'dir')
  mkdir('./bin');
end

% add appropriate folders to path

addpath('./bin');
addpath('./src');
addpath('./examples');

% compile mex binaries

buildqrmex

% test mex binaries

testqrmex

