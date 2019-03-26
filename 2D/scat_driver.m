% driver code for 2D Helmholtz exterior scattering from smooth domain.
% using MFS and FLAM FDS.
% Barnett 3/25/19

clear;
v = 1;        % verbosity (0,1)

k = 30;    % wavenumber

b = 1; a = .3; w = 5; % try 25;        % smooth wobbly radial shape params
R = @(t) b*(1 + a*cos(w*t));

ti = -pi/6; ui = @(x) exp(1i*k*(cos(ti)*real(x)+sin(ti)*imag(x))); % u_inc
x0 = -1 + 1.5i;   % scatt wave soln test pt.  Note pts are in complex notation

N = 1e3;

% MFS: t = surface pt struct, s = source pt struct
M = round(1.2*N);         % # bdry pts
t.t = (1:M)'/M*2*pi; t.x = exp(1i*t.t).*R(t.t);  % bdry pts
s.t = (1:N)'/N*2*pi; s.x = exp(1i*s.t).*R(s.t);  % MFS src pts
imagd = 100/N;     % scale dist w/ 1/N
s.t = 1i*imagd + (1:N)'/N*2*pi; s.x = exp(1i*s.t).*R(s.t);
if v, figure(1); plot(t.x, 'b.-'); hold on; plot(s.x, 'r.'); plot(x0,'+');
  axis equal; title(sprintf('N=%d, imagd=%g',N,imagd)); hold off; drawnow; end

% convert to real coords: rx = surface pts (rows), cx = source pts (cols)...
rx = [real(t.x(:))'; imag(t.x(:))']; cx = [real(s.x(:))'; imag(s.x(:))'];

flampar.rank_or_tol = 1e-12;
p = 64; flampar.p=p; theta = (1:p)*2*pi/p;         % proxy pts vs unit box
flampar.proxy = 1.5*[cos(theta); sin(theta)];
clear theta
flampar.opts = []; flampar.opts.verb = 0;
flampar.occ=128;    % max pts per box for quadtree

lsqpar.tau = eps^(-1/3);  % params for LSQ
lsqpar.qr = 'm';   % 'm'=matlab, 'q'=qr,  's'=spqr (needs SuiteSparse)
lsqpar.refine = 1;      % 0 or 1

meth = 'l';  % 'l'=dense LU, 'q'=dense QR, 'r'=FLAM rskel
F = factor(k,rx,cx,meth,flampar,lsqpar);  % direct solve, into struct fac
rhs = -ui(t.x);     % eval uinc on bdry
co = solve(F,rhs,meth);
nrm = norm(co);
rrms = norm(F.A*co - rhs)/sqrt(M);
%fprintf('rms err in resid = %.3g\n',rrms)

u = mfseval(k,[real(x0);imag(x0)],cx,co,'d');

ue =  -0.764438165507287 +     0.357188677072235i;  % true val (a=.3,w=5,k=30)
disp('      N       resid rms      u err at pt   soln 2-norm:')
Ns = N;
format short g; disp([Ns(:), rrms, abs(u-ue), nrm]); format long g

m = 1e3;    % check at new bdry pts
clear b; b.t = (1:m)'/m*2*pi; b.x = exp(1i*b.t).*R(b.t);  % bdry pts
ub = mfseval(k,[real(b.x)';imag(b.x)'],cx,co,'d');
utotb = ub + ui(b.x).';                 % physical potential on bdry
fprintf('bdry rms err: %.3g \n', norm(utotb)/sqrt(m)) % u_tot on bdry = 0?

if v, dx = 0.03; xmax = 2.0; g = -xmax:dx:xmax;
  [xx yy] = meshgrid(g,g);
  ii = ~insideradfunc(R,xx+1i*yy);
  xx = xx(ii)'; yy = yy(ii)';       % row vecs of ext targ coords only
  tic; uu =  mfseval(k,[xx;yy],cx,co,'d'); fprintf('grid eval %.3g s\n',toc);
  uu = uu + ui([xx + 1i*yy]);   % add incident to get physical potential (rows)
  ug = nan*ii; ug(ii) = uu;     % overwrite ext only
  figure(2); imagesc(g,g,real(ug)); caxis(2*[-1 1]); axis equal tight xy
  hold on; plot(t.x,'-');
end
