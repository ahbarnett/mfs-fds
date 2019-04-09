% driver for 2D Helmholtz exterior scattering from domains w/ 0 or 1 corners,
% and panel discr, using MFS and FLAM FDS.
% Barnett 3/26/19

clear;
v = 1;        % verbosity: 0 = text only, 1 = total potential plot

k = 15;    % wavenumber

dom='c';    % choose one of 3 domain types...
if dom=='s'  % set up Z, Z', Z'' parametrization of chosen domain in [-pi,pi]..
  a = 0.3; w = 3; b = 0.5; t.Z = @(t) (1+a*cos(w*t+b)).*exp(1i*t); % trefoil
elseif dom=='c'   % Bremer teardrop (corner at t=0 for relative accuracy)
  b = pi/2; %2*pi/3; % 3*pi/5; % corner angle: nonreentrant, ie <pi
  a = 2/tan(b/2); t.Z = @(t) -a*sin(abs(t)/2)+1i*sin(t);
elseif dom=='r'  % Bremer boomerang 2pi-b (corner at t=0 for relative accuracy)
  b = pi/2; %<pi; pi/4 needs Nb=20, Nr=90; pi/2 needs Nb=12, Nr=60
  a = 2/tan((2*pi-b)/2)/3; t.Z = @(t) -a*sin(3/2*abs(t))+1i*sin(t);
end

ti = -5*pi/6; ui = @(x) exp(1i*k*(cos(ti)*real(x)+sin(ti)*imag(x))); % u_inc
x0 = -1 + 1.5i;   % scatt wave soln test pt.  Note pts are in complex notation

evalmeth = 'f';   % summation method for pot eval: 'd' direct slow, 'f' FMM.

p = 12;       % src pts per panel
Nb = 20; Nr=40; % defaults: # base panels (Np= # panels), # levels of r-adic ref
if dom=='s', Nr=0; end
reffac = 3;   % refinement factor to corners
o.qtype = 'u';
[t ts pans] = panelquadr(t,p,Nb,Nr,reffac,0,o); N = numel(t.x); Np = numel(ts);

pmfs = p;   % # MFS src pts per pan
o.reldist = -0.3;   % interior sources
s = panelmfs(p,t,ts,pans,o); M = numel(s.x);  % MFS src pts (s = struct)

if v, figure(1); clf; plot(t.x, 'b.-'); hold on; plot(s.x, 'r.'); plot(x0,'+');
  axis equal; title(sprintf('N=%d',N)); hold off; drawnow; end

% convert to real coords: rx = surface pts (rows), cx = source pts (cols)...
rx = [real(t.x(:))'; imag(t.x(:))']; cx = [real(s.x(:))'; imag(s.x(:))'];

% linear solver choice...
meth = 'l';  % 'l'=dense LU, 'q'=dense QR (both O(N^3)); 'r'=FLAM rskel

flampar.rank_or_tol = 1e-12;
flampar.p=64;
flampar.opts = []; flampar.opts.verb = 0;
flampar.occ=512;    % max pts per box for quadtree (might need to tune)

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
fprintf('rms err in resid = %.3g\n',rrms)

u = mfseval(k,[real(x0);imag(x0)],cx,co,'d');   % one targ: always use direct
ue =nan;  % dummy

disp('      N       resid rms      u err at pt   soln 2-norm:')
Ns = N;
format short g; disp([Ns(:), rrms, abs(u-ue), nrm]); format long g

pchk = 7;   % check at new bdry pts, this many (neq p) per panel
b = panelquadr(t,pchk,Nb,Nr,reffac,0,o); m = numel(b.x);
ub = mfseval(k,[real(b.x)';imag(b.x)'],cx,co,evalmeth);
utotb = ub + ui(b.x).';                 % physical potential on bdry
if v, figure(3); clf; plot(1:m, real(utotb), '.-'); xlabel('test pt ind');
  ylabel('total pot on bdry'); end
fprintf('bdry rms err: %.3g \n', norm(utotb)/sqrt(m)) % u_tot on bdry = 0?

if v, dx = 0.01; gx = -4:dx:1; gy = -2:dx:2;    % image plot grid
  [xx yy] = meshgrid(gx,gy);
  ii = true(size(xx)); % ~insideradfunc(R,xx+1i*yy);
  xx = xx(ii)'; yy = yy(ii)';       % row vecs of ext targ coords only
  tic; uu =  mfseval(k,[xx;yy],cx,co,evalmeth);
  fprintf('grid eval %.3g s\n',toc);
  uu = uu + ui([xx + 1i*yy]);   % add incident to get physical potential (rows)
  ug = nan*ii; ug(ii) = uu;     % overwrite ext only
  figure(2); clf;
  imagesc(gx,gy,real(ug)); caxis(2*[-1 1]); axis equal tight xy
  colorbar; hold on; plot(t.x,'-');
end
