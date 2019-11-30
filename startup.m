% matlab set-up paths etc for mfs-fds.
% Barnett 4/8/19

% Path setup --------------------------------------------------------
h = fileparts(mfilename('fullpath'));        % direc of this file
addpath(h);
addpath([h '/fds/common'])
addpath([h '/fds/rskelfm'])
addpath([h '/utils']);
if 0
  addpath(genpath(h))                        % gives access to all subdirs
  rmpath(genpath(fullfile(h,'.git')))
end

addpath([h '/../fmmlib2d/matlab']);            % 2D FMM

% following is correct if FLAM placed alongside mfs-fds directory...
FLAM = [h '/../FLAM'];      % user: adjust to top of your FLAM installation
run([FLAM '/startup.m']);

% following is correct if SuiteSparse placed alongside mfs-fds directory
addpath('../SuiteSparse/SPQR/MATLAB')



% Other setup ------------------------------------------------------------

maxNumCompThreads('automatic');   % gives 1 or 2 thr/core, depending on machine
%maxNumCompThreads(2*maxNumCompThreads('automatic'));  % bump to 2 thr/core?
fprintf('MATLAB will use %d threads\n',maxNumCompThreads)

set(0,'showHiddenHandles','on');        % eg lines linesmoothing

format long g
format compact

set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
                    .75 0 .75; .75 .75 0; .25 .25 .25])
set(groot,'DefaultFigureColormap',jet(256))

