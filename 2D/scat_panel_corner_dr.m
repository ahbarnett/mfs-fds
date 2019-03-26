% driver code for 2D Helmholtz exterior scattering from domains w/ 0,1 corner,
% and panel discr, using MFS and FLAM FDS.
% Barnett 3/26/19

clear;
v = 1;        % verbosity (0,1)
k = 15;    % wavenumber
dom='c';
if dom=='s'  % set up Z, Z', Z'' parametrization of chosen domain in [-pi,pi]..
  a = 0.3; w = 3; b = 0.5; t.Z = @(t) (1+a*cos(w*t+b)).*exp(1i*t); % trefoil
elseif dom=='c'   % Bremer teardrop (corner at t=0 for relative accuracy)
  b = pi/3; %2*pi/3; % 3*pi/5; % corner angle: nonreentrant, ie <pi
  a = 2/tan(b/2); t.Z = @(t) -a*sin(abs(t)/2)+1i*sin(t);
elseif dom=='r'  % Bremer boomerang 2pi-b (corner at t=0 for relative accuracy)
  b = pi/2; %<pi; pi/4 needs Nb=20, Nr=90; pi/2 needs Nb=12, Nr=60
  a = 2/tan((2*pi-b)/2)/3; t.Z = @(t) -a*sin(3/2*abs(t))+1i*sin(t);
end

ti = -5*pi/6; ui = @(x) exp(1i*k*(cos(ti)*real(x)+sin(ti)*imag(x))); % u_inc
x0 = -1 + 1.5i;   % scatt wave soln test pt.  Note pts are in complex notation

p = 12;   % pts per panel

Nb = 20; Nr=20; % defaults: # base panels (Np= # panels), # levels of r-adic ref
if dom=='s', Nr=0; end
reffac = 3;   % refinement factor to corners
o.qtype = 'u';
[t ts pans] = panelquadr(t,p,Nb,Nr,reffac,0,o); N = numel(t.x); Np = numel(ts);

pmfs = p;   % # MFS pts per pan
o.reldist = -0.3;   % interior sources
s = panelmfs(p,t,ts,pans,o); M = numel(s.x);  % MFS src pts (s = struct)

if v, figure(1); clf; plot(t.x, 'b.-'); hold on; plot(s.x, 'r.'); plot(x0,'+');
  axis equal; title(sprintf('N=%d',N)); hold off; drawnow; end

% convert to real coords: rx = surface pts (rows), cx = source pts (cols)...
rx = [real(t.x(:))'; imag(t.x(:))']; cx = [real(s.x(:))'; imag(s.x(:))'];

flampar.rank_or_tol = 1e-12;
flampar.p=64;
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
fprintf('rms err in resid = %.3g\n',rrms)

u = mfseval(k,[real(x0);imag(x0)],cx,co,'d');
ue =nan;

disp('      N       resid rms      u err at pt   soln 2-norm:')
Ns = N;
format short g; disp([Ns(:), rrms, abs(u-ue), nrm]); format long g

pchk = 7;   % check at new bdry pts, this many per panel
b = panelquadr(t,pchk,Nb,Nr,reffac,0,o); m = numel(b.x);
ub = mfseval(k,[real(b.x)';imag(b.x)'],cx,co,'d');
utotb = ub + ui(b.x).';                 % physical potential on bdry
fprintf('bdry rms err: %.3g \n', norm(utotb)/sqrt(m)) % u_tot on bdry = 0?

if v, dx = 0.03; gx = -4:dx:1; gy = -2:dx:2;
  [xx yy] = meshgrid(gx,gy);
  ii = true(size(xx)); % ~insideradfunc(R,xx+1i*yy);
  xx = xx(ii)'; yy = yy(ii)';       % row vecs of ext targ coords only
  tic; uu =  mfseval(k,[xx;yy],cx,co,'d'); fprintf('grid eval %.3g s\n',toc);
  uu = uu + ui([xx + 1i*yy]);   % add incident to get physical potential (rows)
  ug = nan*ii; ug(ii) = uu;     % overwrite ext only
  figure(2); clf;
  imagesc(gx,gy,real(ug)); caxis(2*[-1 1]); axis equal tight xy
  hold on; plot(t.x,'-');
end
