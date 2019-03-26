% matlab set-up for mfs-fds. Barnett 3/25/19

dirs = {'.','2D'};
for s = dirs, addpath(sprintf('%s',s{:})); end, clear s

FLAM = '../FLAM';            % user adjust to top of FLAM installation
run([FLAM '/startup.m']);

