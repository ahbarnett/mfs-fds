% driver for 2D Helmholtz exterior scattering from domains w/ 0 or 1 corners,
% and panel discr, using MFS and FLAM FDS, with convergence & accurate answers.
% Barnett 3/26/19; "exact" values to compare far-field 4/22/19, 5 domains.

v = 1;        % verbosity: 0 = text only, 1 = total potential plot

k = 15;    % wavenumber (uex for k=15 and the below x0 & uinc)

o.reldist = -0.5; o.meth = 'c';    % default MFS interior source loc algorithm

if ~exist('dom','var'), dom='2'; end      % choose one of 5 domain types...
if dom=='s'  % set up Z, Z', Z'' parametrization of chosen domain in [-pi,pi]..
  a = 0.3; w = 3; b = 0.5; t.Z = @(t) (1+a*cos(w*t+b)).*exp(1i*t); % trefoil
  Nb = 26; Nr=0;   %  # base panels (Np= # panels), # levels of r-adic ref
  reffac = nan;
  uex = 0.0115203463698269 -     0.802730582580788i;  % to maybe 1e-14
elseif dom=='2'   % Bremer teardrop (corner at t=0 for relative accuracy)
  b = pi/2;     % corner angle: reentrant, ie >pi, for ext BVP
  a = 2/tan(b/2); t.Z = @(t) -2*sin(abs(t)/2)+(2i/a)*sin(t);
  Nb = 30; Nr=40;   %  # base panels (Np= # panels), # levels of r-adic ref
  reffac = 3;   % refinement factor to corners
  uex = 0.0425178324367369 -      0.19364796049328i;   % to maybe 1e-14
elseif dom=='4'   % Bremer teardrop,  narrower angle (still length=2)
  b = pi/4;       % corner angle: reentrant, ie >pi, for ext BVP
  a = 2/tan(b/2); t.Z = @(t) -2*sin(abs(t)/2)+(2i/a)*sin(t);
  Nb = 40; Nr=40;   %  # base panels (Np= # panels), # levels of r-adic ref
  reffac = 2.5;   % refinement factor to corners
  o.reldist = -0.45;   % override
  uex = 0.0320438849594651 -     0.103076591191477i;    % to maybe 1e-14
elseif dom=='8'   % Bremer teardrop,  even narrower angle (still length=2)
  b = pi/8;       % corner angle: reentrant, ie >pi, for ext BVP
  a = 2/tan(b/2); t.Z = @(t) -2*sin(abs(t)/2)+(2i/a)*sin(t);
  Nb = 40; Nr=50;   %  # base panels (Np= # panels), # levels of r-adic ref
  reffac = 1.8;   % refinement factor to corners
  o.reldist = -0.3;   % override
  uex = 0.0282227868146729 -    0.0889087732747555i;    % to maybe 1e-11
elseif dom=='r'  % Bremer boomerang 2pi-b (corner at t=0 for relative accuracy)
  b = pi/2;    % nonreentrant >pi, not a challenge
  a = 2/tan((2*pi-b)/2)/3; t.Z = @(t) -a*sin(3/2*abs(t))+1i*sin(t);
  Nb = 30; Nr=10;   %  # base panels (Np= # panels), # levels of r-adic ref
  reffac = 3;   % refinement factor to corners
  uex = 0.506135119824532 +       1.3602391736054i;      % to maybe 1e-14
end

ti = -5*pi/6; ui = @(x) exp(1i*k*(cos(ti)*real(x)+sin(ti)*imag(x))); % u_inc
x0 = 1.5 + 0.9i;   % scatt wave soln test pt (complex notation); don't change

evalmeth = 'f';   % summation method for pot eval: 'd' direct slow, 'f' FMM.

pmfs = 14;    % # MFS src pts per pan
p = round(1.2*pmfs); o.qtype = 'u';      % surf pts per panel

mode = 'B';   % conv mode: '1' single run; see below for others...
Nbs = Nb; Nrs = Nr; c = (-5:0);   % conv up to the given value, will get scaled
if mode=='N', Nbs = Nb+2*c; Nrs = Nr + 0*Nbs;    % Nb conv only
elseif mode=='C', Nrs = max(0,Nr+5*c); Nbs = Nb + 0*Nrs; % Nr conv only
elseif mode=='B', Nbs = Nb+2*c; Nrs = max(0,Nr+5*c); end   % both conv together

us = nan*Nbs; rrms=us; nrm=us; brms=us;  % things to track the conv of
for i=1:numel(Nbs),  Nb=Nbs(i); Nr=Nrs(i); % ----------------- conv loop ----

  [t tpan] = panelquadr_1corn(t,p,Nb,Nr,reffac,0,o);
  M = numel(t.x); Np = numel(tpan);

  s = panelcurvemfs(pmfs,t,tpan,o); N = numel(s.x);  % MFS src pts (s = struct)
  
  if v, figure(1); clf; plot(t.x, 'b.-'); hold on; plot(s.x, 'r.');
    plot(x0,'+'); axis equal; title(sprintf('N=%d',N)); hold off; drawnow;
  end

  % convert to real coords: rx = surface pts (rows), cx = source pts (cols)...
  rx = [real(t.x(:))'; imag(t.x(:))']; cx = [real(s.x(:))'; imag(s.x(:))'];

  % linear solver choice... (these are small so dense is ok)
  meth = 'l';  % 'l'=dense LU, 'q'=dense QR (both O(N^3)); 'r'=FLAM rskel
  
  flampar.rank_or_tol = 1e-12;
  flampar.p=64;
  flampar.opts = []; flampar.opts.verb = 0;
  flampar.occ=512;    % max pts per box for quadtree (might need to tune)

  lsqpar.tau = eps^(-1/3);  % params for LSQ
  lsqpar.qr = 'q';   % 'q'=qr, 's'=spqr (needs SuiteSparse)
  lsqpar.refine = 0;      % 0 or 1

  F = factor(k,rx,cx,meth,flampar,lsqpar);  % direct solve, into struct fac
  rhs = -ui(t.x);     % eval uinc on bdry
  co = solve(F,rhs,meth);
  nrm(i) = norm(co);
  ur = mfseval(k,[real(t.x)';imag(t.x)'],cx,co,evalmeth);  % apply A to co
  rrms(i) = norm(ur(:) - rhs)/sqrt(M);

  us(i) = mfseval(k,[real(x0);imag(x0)],cx,co,'d');  % 1 targ: always use direct

  fprintf('pmfs=%d, p=%d, Nb=%d, Nr=%d: resid rms = %.3g, co 2-nrm = %.3g, Re(u) = %.15g\n',pmfs,p,Nb,Nr,rrms(i),nrm(i),real(us(i)))

  pchk = 7;   % check at new bdry pts, this many (neq p) per panel
  b = panelquadr_1corn(t,pchk,Nb,Nr,reffac,0,o); m = numel(b.x);
  ub = mfseval(k,[real(b.x)';imag(b.x)'],cx,co,evalmeth);
  utotb = ub + ui(b.x).';                 % physical potential on bdry
  brms(i) = norm(utotb)/sqrt(m);
  fprintf('\tbdry rms err: %.3g \t est u err at pt = %.3g\n',brms(i),abs(us(i)-uex)) % u_tot on bdry = 0?

end                      % -------------------- end conv ----------------

% show conv...
fprintf('Nb\tNr    resid rms  co 2-nrm  bdry rms  est u err  Re(u) :\n');
for i=1:numel(Nbs), fprintf('%d\t%d\t%.2g\t%.2g\t%.2g\t%.3g   \t%.16g\n', Nbs(i),Nrs(i),rrms(i), nrm(i),brms(i),abs(us(i)-uex),real(us(i))); end
figure(4); clf; Nsho=Nbs; ssho='N_b';
if mode=='C', Nsho=Nrs; ssho='N_r'; elseif mode=='B',ssho='N_b (& N_r incr)';end
semilogy(Nsho,[abs(us-uex);rrms;brms],'+-'); xlabel(ssho);
ylabel('err'); legend('u far pt est err','resid rms','bdry rms');
title('soln conv');

if v, figure(3); clf; subplot(2,1,1); semilogy(1:m, abs(utotb), '.-'); % err
  xlabel('test pt ind'); ylabel('total pot on bdry'); axis tight;
  subplot(2,1,2); semilogy(1:N, abs(co), '.-'); xlabel('src pt ind j');
  ylabel('coeff_j'); axis tight;
end

if v, dx = 0.01; gx = -3:dx:2; gy = -2:dx:2;    % image plot grid, covers all 3
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
