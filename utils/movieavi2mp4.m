function status = movieavi2mp4(nam,fps)
% MOVEAVI2MP4.  Use system call to ffmpeg to convert AVI to MP4
%
% movieavi2mp4('f') converts f.avi to f.mp4, overwriting the latter.
%
% movieavi2mp4('f',30) does the same, fixing 30 FPS in the MP4.
%
% Returns the status. Linux systems only.
% There is some code to customize to various of our linux environments.

% Barnett 4/19/19

if nargin<2, fps = 20; end   % default frames per sec

[~,kernel] = system('uname -r');
cmd = sprintf('ffmpeg -i %s.avi -y -c:v libx264 -crf 20 -r %d %s.mp4',...
              nam,fps,nam);

if strcmp(kernel,'el7')            % on flatiron cluster
  disp('running on EL7 system - will module load ffmpeg...')
  cmd = ['module load lib/x264; module load ffmpeg; ' cmd];
else                               % on generic linux
  disp('running on generic system - will remove matlab ld_path...')
  cmd = ['unset LD_LIBRARY_PATH; ' cmd];
end

status = system(cmd);  % do it
