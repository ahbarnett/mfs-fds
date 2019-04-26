% driver code for 2D Helmholtz exterior scattering from smooth domain.
% using MFS and FLAM FDS. Small example to test: FLAM and FMM2D.

% Barnett 3/25/19. 4/8/19: added FMM2D

clear;
v = 1;        % verbosity: 0 = text only, 1 = total potential plot

k = 30;    % wavenumber (fix all params if want to compare to ue value below)

b = 1; a = .3; w = 5; % try 25;        % smooth wobbly radial shape params
R = @(t) b*(1 + a*cos(w*t));

ti = -pi/6; ui = @(x) exp(1i*k*(cos(ti)*real(x)+sin(ti)*imag(x))); % u_inc
x0 = -1 + 1.5i;   % scatt wave soln test pt.  Note pts are in complex notation

N = 1e4;         % overkill by factor 20, but good for test evalmeth='f'.
evalmeth = 'f';   % summation method for pot eval: 'd' direct slow, 'f' FMM.

% MFS: t = surface pt struct, s = source pt struct
M = round(1.2*N);         % # bdry pts
t.t = (1:M)'/M*2*pi; t.x = exp(1i*t.t).*R(t.t);  % bdry pts
s.t = (1:N)'/N*2*pi; s.x = exp(1i*s.t).*R(s.t);  % MFS src pts
imagd = 0.1;                   % Note for imagd=0.1 in this BVP, rank(A)<600
s.t = 1i*imagd + (1:N)'/N*2*pi; s.x = exp(1i*s.t).*R(s.t);
fprintf('min src ppw = %.3g; bdry ppw = %.3g\n',2*pi/(k*max(abs(diff(s.x)))),2*pi/(k*max(abs(diff(t.x)))))
if v, figure(1); clf; plot(t.x, 'b.-'); hold on; plot(s.x, 'r.'); plot(x0,'+');
  axis equal; title(sprintf('N=%d, imagd=%g',N,imagd)); hold off; drawnow; end

% convert to real coords: rx = surface pts (rows), cx = source pts (cols)...
rx = [real(t.x(:))'; imag(t.x(:))']; cx = [real(s.x(:))'; imag(s.x(:))'];

% linear solver choice...
meth = 'r';  % 'l'=dense LU, 'q'=dense QR (both O(N^3)); 'r'=FLAM rskel

flampar.rank_or_tol = 1e-12;
flampar.p=64;  % num proxy pts (should depend on eps)
flampar.opts = []; flampar.opts.verb = 0;
flampar.occ=128;    % max pts per box for quadtree

lsqpar.tau = eps^(-1/3);  % params for LSQ
lsqpar.qr = 'q';   % 'q'=qr, 's'=spqr (needs SuiteSparse)
lsqpar.refine = 0;      % 0 or 1

F = factor(k,rx,cx,meth,flampar,lsqpar);  % direct solve, into struct fac
if meth~='r', fprintf('eps-rank of A : %d\n',rank(F.A,1e-14*normest(F.A))), end
rhs = -ui(t.x);     % eval uinc on bdry
co = solve(F,rhs,meth);
nrm = norm(co);
ur = mfseval(k,[real(t.x)';imag(t.x)'],cx,co,evalmeth);  % apply A to co
rrms = norm(ur(:) - rhs)/sqrt(M);

u = mfseval(k,[real(x0);imag(x0)],cx,co,'d');   % one targ: always use direct

ue =  -0.764438165507287 +     0.357188677072235i;  % true val (a=.3,w=5,k=30)
disp('      N       resid rms      u err at pt   soln 2-norm:')
Ns = N;
format short g; disp([Ns(:), rrms, abs(u-ue), nrm]); format long g

m = 1e3;    % check at this many new shifted-grid bdry pts
clear b; b.t = (0.5:m-0.5)'/m*2*pi; b.x = exp(1i*b.t).*R(b.t);
ub = mfseval(k,[real(b.x)';imag(b.x)'],cx,co,evalmeth);
utotb = ub + ui(b.x).';                 % physical potential on bdry
fprintf('bdry rms err: %.3g \n', norm(utotb)/sqrt(m)) % u_tot on bdry = 0?

if v, dx = 0.01; xmax = 2.0; g = -xmax:dx:xmax;      % image plot grid
  [xx yy] = meshgrid(g,g);
  ii = ~insideradfunc(R,xx+1i*yy);
  xx = xx(ii)'; yy = yy(ii)';       % row vecs of ext targ coords only
  tic; uu = mfseval(k,[xx;yy],cx,co,evalmeth);
  fprintf('grid eval %.3g s\n',toc);
  uu = uu + ui([xx + 1i*yy]);   % add incident to get physical potential (rows)
  ug = nan*ii; ug(ii) = uu;     % overwrite ext only
  figure(2); clf; imagesc(g,g,real(ug)); caxis(2*[-1 1]); axis equal tight xy
  colorbar; hold on; plot(t.x,'-');
end
