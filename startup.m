% matlab set-up paths etc for mfs-fds.
% Barnett 3/26/19

h = fileparts(mfilename('fullpath'));        % direc of this file
addpath(h);
if 0
  addpath(genpath(h))                        % gives access to all subdirs
  rmpath(genpath(fullfile(h,'.git')))
end
  
% correct if FLAM placed alongside mfs-fds directory...
FLAM = [h '/../FLAM'];            % user: adjust to top of FLAM installation
run([FLAM '/startup.m']);

% double # threads (hyperthreading): accels some multithreaded matlab ops!...
%maxNumCompThreads(2*maxNumCompThreads('automatic'));
maxNumCompThreads('automatic');
fprintf('MATLAB will use %d threads\n',maxNumCompThreads)

set(0,'showHiddenHandles','on');        % eg lines linesmoothing

format long g
format compact

set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
                    .75 0 .75; .75 .75 0; .25 .25 .25])
set(groot,'DefaultFigureColormap',jet(256))

% SuiteSparse paths
SuiteSparse = '../SuiteSparse';
addpath([SuiteSparse])
addpath([SuiteSparse '/UMFPACK/MATLAB'])
addpath([SuiteSparse '/CHOLMOD/MATLAB'])
addpath([SuiteSparse '/AMD/MATLAB'])
addpath([SuiteSparse '/COLAMD/MATLAB'])
addpath([SuiteSparse '/CCOLAMD/MATLAB'])
addpath([SuiteSparse '/CAMD/MATLAB'])
addpath([SuiteSparse '/ssget'])
addpath([SuiteSparse '/CXSparse/MATLAB/Demo'])
addpath([SuiteSparse '/CXSparse/MATLAB/CSparse'])
addpath([SuiteSparse '/LDL/MATLAB'])
addpath([SuiteSparse '/BTF/MATLAB'])
addpath([SuiteSparse '/KLU/MATLAB'])
addpath([SuiteSparse '/SPQR/MATLAB'])
addpath([SuiteSparse '/RBio/RBio'])
addpath([SuiteSparse '/MATLAB_Tools'])
addpath([SuiteSparse '/MATLAB_Tools/Factorize'])
addpath([SuiteSparse '/MATLAB_Tools/MESHND'])
addpath([SuiteSparse '/MATLAB_Tools/LINFACTOR'])
%addpath([SuiteSparse '/MATLAB_Tools/findcomponents'])
addpath([SuiteSparse '/MATLAB_Tools/GEE'])
addpath([SuiteSparse '/MATLAB_Tools/shellgui'])
addpath([SuiteSparse '/MATLAB_Tools/waitmex'])
addpath([SuiteSparse '/MATLAB_Tools/spqr_rank'])
addpath([SuiteSparse '/MATLAB_Tools/spqr_rank/SJget'])
addpath([SuiteSparse '/MATLAB_Tools/SuiteSparseCollection'])
addpath([SuiteSparse '/MATLAB_Tools/SSMULT'])
addpath([SuiteSparse '/MATLAB_Tools/dimacs10'])
addpath([SuiteSparse '/MATLAB_Tools/spok'])
addpath([SuiteSparse '/MATLAB_Tools/sparseinv'])
addpath([SuiteSparse '/Mongoose/MATLAB'])
