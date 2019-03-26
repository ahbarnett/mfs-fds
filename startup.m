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
maxNumCompThreads(2*maxNumCompThreads('automatic'));

set(0,'showHiddenHandles','on');        % eg lines linesmoothing

format long g
format compact

set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
                    .75 0 .75; .75 .75 0; .25 .25 .25])
set(groot,'DefaultFigureColormap',jet(256))
