% driver code for 3D Helmholtz exterior scattering from smooth domain.
% using MFS and FLAM FDS. Small example to test: FLAM and FMM3D.

% Liu 2/27/21

clear;
tic; 
filename = 'FDS_3D_ellipsoid.txt';
fid = fopen(filename, 'a');
k = 3;    % wavenumber (fix all params if want to compare to ue value below)
theta = -pi/4;  phi = pi/3; 
kx = k * cos(theta) * cos(phi);  ky = k * cos(theta) * sin(phi);  kz = k * sin(theta); 

N = 1e4;         % overkill by factor 20, but good for test evalmeth='f'.
evalmeth = 'f';   % summation method for pot eval: 'd' direct slow, 'f' FMM
M = 1.2 * N; 
P = 100; 
N0 = N/P; 
M0 = M/P; 
% MFS: t = surface pt struct, s = source pt struct
d = 0.1;                   % Note for imagd=0.1 in this BVP, rank(A)<600
Rx = 0.6; Ry = 0.6; Rz = 1.0; 
t = getShapeEllip(Rx, Ry, Rz, M0, P);
s = getSourceEllip(Rx, Ry, Rz, d, N0, P);
f = QuasiPerVecFilling3D(t, kx, ky, kz);

rx = [t.x'; t.y'; t.z']; cx = [s.x'; s.y'; s.z'];

% linear solver choice...
meth = 'r';  % 'l'=dense LU, 'q'=dense QR (both O(N^3)); 'r'=FLAM rskel

flampar.rank_or_tol = 1e-10;
flampar.rp = 1.5; % radius of proxy sphere 
flampar.p = 32;   % num proxy pts (should depend on eps)
flampar.opts = []; flampar.opts.verb = 0;
flampar.occ = 128;    % max pts per box for quadtree

fprintf('params: N = %d, Rx = %.3g, Ry = %.3g, Rz = %.3g, k = %d, d = %.3g.\n',N, Rx, Ry, Rz, k, d); 
fprintf(fid, 'params: N = %d, Rx = %.3g, Ry = %.3g, Rz = %.3g, k = %d, d = %.3g.\n',N, Rx, Ry, Rz, k, d); 

fprintf('FLAM params: proxy R = %.3g, tol = %e, p = %d, occ = %d.\n', flampar.rp, flampar.rank_or_tol, flampar.p, flampar.occ); 
fprintf(fid, 'FLAM params: proxy R = %.3g, tol = %e, p = %d, occ = %d.\n', flampar.rp, flampar.rank_or_tol, flampar.p, flampar.occ); 

% params for LSQ
lsqpar.meth = 'u';  % 'u'=underdetermined, 'o'=overdetermined
lsqpar.tau = eps^(-1/3);  % constraint weighting
lsqpar.lambda = 1e-6;  % regularization in OLS
lsqpar.qr = 'q';   % 'q'=qr, 's'=spqr (needs SuiteSparse)
lsqpar.refine = 0;      % 0 or 1

F = factor(k,rx,cx,meth,flampar,lsqpar);  % direct solve, into struct fac
if meth~='r', fprintf('eps-rank of A : %d\n',rank(F.A,1e-14*normest(F.A))), end
rhs = f;     % eval uinc on bdry
co = solve(F,rhs,meth);
nrm = norm(co);
ur = mfseval(k,[t.x'; t.y'; t.z'],cx,co,evalmeth);  % apply A to co
rrms = norm(ur(:) - rhs)/sqrt(M);

m = 37;    % check at this many new shifted-grid bdry pts
b = getShapeEllip(Rx, Ry, Rz, m, m);
ub = mfseval(k,[b.x'; b.y'; b.z'],cx,co,evalmeth);
uib = QuasiPerVecFilling3D(b, kx, ky, kz);
utotb = ub(:) - uib;                 % physical potential on bdry
fprintf('bdry rms err: %.3g \n', norm(utotb)/m); % u_tot on bdry = 0?
fprintf(fid, 'bdry rms err: %.3g \n', norm(utotb)/m); % u_tot on bdry = 0?
fprintf(fid, 'the whold program takes %.3g seconds.\n', toc);
toc; 
fprintf(fid, '------------------------------------------------------------ \n\n');